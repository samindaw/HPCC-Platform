/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include <cstdlib>
#include <queue>
#include <typeinfo>
#include "mpi_wrapper.hpp"
#include "mputil.hpp"

#define WAIT_DELAY 100
#define PROBING_DELAY 100

#define DUAL_MESSAGE_SCHEME
#define USE_BLOCK_SEND
//#define USE_BLOCK_RECV

//----------Functions and Data structures managing communication data related to Send/Recv Communications in orogress-----------//

class ControlMessage
{
public:
    int dataTag;
    mptag_t replyTag;
    ControlMessage(){}
    ControlMessage(int _dataTag,mptag_t _replyTag): dataTag(_dataTag), replyTag(_replyTag){}
};

// Data structure to keep the data relating to send/receive communications
class CommData
{
private:
    bool cancellationLock = false;
    bool cancellationInProgress = false;
    bool completed = false;
    MPI_Request req;           // persistent request object to keep track of ongoing MPI call
    const hpcc_mpi::MPIComm& comm;                // MPI communicator

    CriticalSection dataChangeLock;
    void lock(){dataChangeLock.enter();}
    void unlock(){dataChangeLock.leave();}

public:    
    int rank;                       // source/destination rank of the processor
    int tag;                        // MPI tag information
    ControlMessage cm;
    bool isEqual(int _rank, int _tag)
    {
        return (_rank==MPI_ANY_SOURCE || rank==_rank) && (_tag==MPI_ANY_TAG || tag==_tag);
    }

    CommData(int _rank, int _tag, const hpcc_mpi::MPIComm& _comm):
        rank(_rank), tag(_tag), comm(_comm){}

    MPI_Request *request()
    {
        return &req;
    }
    bool notifyCancellation()
    {
        CriticalBlock block(dataChangeLock);
        if (!completed)
            cancellationInProgress = true;
        return !completed;
    }
    bool isCompleted()
    {
        CriticalBlock block(dataChangeLock);
        return completed;
    }
    void setCompleted()
    {
        CriticalBlock block(dataChangeLock);
        completed = true;
    }

    bool lockFromCancellation()         //returns true only if currently some thread is not attempting to cancel
    {
        lock();
        if (cancellationInProgress)
        {
            unlock();
        } else
        {
            cancellationLock = true;
        }
        return !cancellationInProgress;
    }

    void releaseCancellationLock()
    {
        if (cancellationLock)
        {
            cancellationLock = false;
            unlock();
        }//TODO else clause: throw a meaningful error for invalid call for this function
    }

    ~CommData(){}
};

std::vector<CommData*> asyncCommData; // CommData list to manage while send/recv communication in progress
CriticalSection commDataLock;         // A mutex lock for the index list above

CommData* _popCommData(size_t index)
{
    CommData* ret = NULL;
    assertex(index>=0);
    assertex(index < asyncCommData.size());
    ret = asyncCommData[index];
    asyncCommData.erase(asyncCommData.begin() + index);
    return ret;
}

void addCommData(CommData *commData)
{
    CriticalBlock block(commDataLock);
    asyncCommData.push_back(commData);
}

void popCommData(int rank, int tag, bool completed, const hpcc_mpi::MPIComm& comm, std::vector<CommData*> &matchedResults)
{
    CriticalBlock block(commDataLock);
   size_t size = asyncCommData.size();
   for(size_t i=(size-1); i>=0; i--)
   {
       if ((asyncCommData[i]->isCompleted()==completed) && (asyncCommData[i]->isEqual(rank, tag)))
           matchedResults.push_back(_popCommData(i));
   }
}

CommData* popCommData(CommData * commData)
{
    CriticalBlock block(commDataLock);
    size_t size = asyncCommData.size();
    for(size_t i=(size-1); i>=0; i--)
    {
        if (asyncCommData[i]==commData)
        {
            _popCommData(i);
            break;
        }
    }
    return commData;
}

//----------------------------------------------------------------------------------------------------------------------------//

int TAG_UB = (1<<28); //INT_MAX //MPI_TAG_UB
int TAG_RANGE = (TAG_UB/2);

//Control message tag values
int TAG_CONTROL_START = 0;
int TAG_CONTROL_UB = TAG_RANGE;
int TAG_CONTROL_SHIFT_BASE = (TAG_CONTROL_UB / 2); // reserve half range for special tags (inc. reply tags)

//Data message tag values
int TAG_DATA_START  = (TAG_CONTROL_UB + 1);
int TAG_DATA_UB  = TAG_UB;
int TAG_DATA_RANGE = (TAG_UB - TAG_RANGE);

CriticalSection dataTagBlock;
int nextAvaialbleDataTag = TAG_DATA_START;
int getNextDataTag()
{
    //Assumption:   We will consume tags in round-robin fasion assuming when a tag
    //              is reused the earlier communication relating to the tag has being completed
    CriticalBlock block(dataTagBlock);
    int next = nextAvaialbleDataTag;
    nextAvaialbleDataTag++;
    if (nextAvaialbleDataTag == TAG_DATA_UB)
        nextAvaialbleDataTag = TAG_DATA_START;
    return next;
}

int getRank(rank_t sourceRank)
{
    if (sourceRank == RANK_ALL)
        return MPI_ANY_SOURCE;
    else
        return sourceRank;
}

int getMPITag(mptag_t mptag)
{

    if (mptag == TAG_ALL)
        return MPI_ANY_TAG;
    else
    {
        int tag = (int)mptag;
        if (tag < 0)
        {
            assertex(-tag < TAG_CONTROL_UB-TAG_CONTROL_SHIFT_BASE); // check not out of range
            return -tag + TAG_CONTROL_SHIFT_BASE;
        }
        else
        {
            assertex(tag <= TAG_CONTROL_UB);
            return tag;
        }
    }
}

mptag_t getMPTag(int tag)
{
    unsigned returnVal;
    if (tag == MPI_ANY_TAG)
        returnVal = TAG_ALL;
    else
    {
        if (tag > TAG_CONTROL_SHIFT_BASE)
            returnVal = (mptag_t) -(tag-TAG_CONTROL_SHIFT_BASE);
        else
            returnVal = tag;
    }
    return (mptag_t)returnVal;
}

MPI_Status waitToComplete( bool& completed, bool& error, bool& canceled, bool& timedout, unsigned timeout, CommData *commData)
{
    CTimeMon tm(timeout);
    MPI_Status stat;
    unsigned remaining;
    bool noCancellation = true;
    int flag;
    while (!(completed || error || (timedout = tm.timedout(&remaining))) && (noCancellation = commData->lockFromCancellation()))
    {

        MPI_Test(commData->request(), &flag, &stat);
        completed = flag > 0;
        commData->releaseCancellationLock();
        usleep(WAIT_DELAY);
    }
    if (completed)
    {
        MPI_Test_cancelled( &stat, &flag);
        canceled = flag > 0;
    }
    else
        canceled = !noCancellation;

    return stat;
}

MPI_Status hasIncomingData(int sourceRank, int mptag, hpcc_mpi::MPIComm& comm,
        bool& incomingMessage, bool& error, bool& canceled, unsigned timeout, CommData *commData)
{

    CTimeMon tm(timeout);
    MPI_Status stat;

    int flag;
    unsigned remaining;
    bool noCancellation = true;
    while (!(incomingMessage || error || (timeout !=0 && tm.timedout(&remaining))) && (noCancellation = commData->lockFromCancellation()))
    {
        MPI_Iprobe(sourceRank, mptag, comm.get(), &flag, &stat);
        incomingMessage = flag > 0;
        commData->releaseCancellationLock();
        if (timeout == 0) break;
        usleep(PROBING_DELAY);
    }
    canceled = !noCancellation;
    return stat;
}

//----------------------------------------------------------------------------//

/** See mpi_wrapper.hpp header file for function descriptions of the following **/

bool hpcc_mpi::MPIComm::hasIncomingMessage(rank_t &sourceRank, mptag_t &mptag)
{
    mpiInitializedCheck();
    bool incomingMessage = false; bool error = false; bool canceled = false;
    int source = getRank(sourceRank);
    int tag = getMPITag(mptag);
    CommData* tmpCommData = new CommData(source, tag, *this);
    MPI_Status stat = hasIncomingData(source, tag, *this, incomingMessage, error, canceled, 0, tmpCommData);
    if (incomingMessage)
    {
        sourceRank = stat.MPI_SOURCE;
        mptag = getMPTag(stat.MPI_TAG);
    }

    return incomingMessage;
}

bool cancelComm(CommData* commData)
{
    bool ret = true;
    if (commData)
    {
        MPI_Status stat;
        int flag;
        //TODO: try catch for error
        MPI_Test(commData->request(),&flag, &stat);
        bool completed = flag > 0;
        if (!completed)
        {
            MPI_Cancel(commData->request());
            MPI_Request_free(commData->request());
        }
    }

    return ret;
}

bool hpcc_mpi::MPIComm::cancelComm(rank_t rank, mptag_t mptag)
{
    mpiInitializedCheck();
    int r = getRank(rank);
    int tag = getMPITag(mptag);
    std::vector<CommData*> commDataList;
    popCommData(r, tag, false, *this, commDataList);
    bool success = true;
    for (auto &commData: commDataList)
    {
        if (commData->notifyCancellation())              // In case readData functions still in progress we want to let
        {                                                                    // them know that we are about to screw up their plans
            if (!::cancelComm(commData))
                success = false;                        // If we managed to cancel everything then the cancellation was successful
            delete commData;
        } // else already recv operation completed
    }
    return success;
}

rank_t hpcc_mpi::MPIComm::rank()
{
    mpiInitializedCheck();
    int rank;
    MPI_Comm_rank(this->get(), &rank);
    return rank;
}

rank_t hpcc_mpi::MPIComm::size()
{
    mpiInitializedCheck();
    int size;
    MPI_Comm_size(this->get(), &size);
    return size;
}

#define MPI_BUFFER_SIZE 16777216        // 16MB
void* mpi_buffer_data;

void hpcc_mpi::initialize()
{
    mpiInitializedCheck(false);

    int required = MPI_THREAD_MULTIPLE;
    int provided;

    MPI_Init_thread(NULL, NULL, required, &provided);
    assertex(provided == required);

    //Allocate memory for the MPI buffered send
    mpi_buffer_data = malloc(MPI_BUFFER_SIZE);
    MPI_Buffer_attach(mpi_buffer_data, MPI_BUFFER_SIZE);
}

void hpcc_mpi::MPIComm::setErrorHandler(MPI_Errhandler handler)
{
    mpiInitializedCheck();
    MPI_Comm_set_errhandler(this->get(),handler);
}

void hpcc_mpi::finalize()
{
    mpiInitializedCheck();
    int size;
    void *ptr;
    MPI_Buffer_detach(&ptr, &size);
    free(mpi_buffer_data);
    MPI_Finalize();
}

void pushMPTag(CMessageBuffer &mbuf, mptag_t tag)
{
    serializeMPtag(mbuf, tag);
}

mptag_t popMPTag(CMessageBuffer &mbuf)
{
    size_t currentPos = mbuf.getPos();
    size_t newLength = mbuf.length() - sizeof(mptag_t);
    mbuf.reset(newLength);
    mptag_t replyTag;
    deserializeMPtag(mbuf, replyTag);
    mbuf.setLength(newLength);
    if (currentPos > newLength)
    {
        currentPos = newLength;
    }
    mbuf.reset(currentPos);
    return replyTag;
}

hpcc_mpi::CommStatus hpcc_mpi::MPIComm::sendData(rank_t dstRank, mptag_t mptag, CMessageBuffer &mbuf, unsigned timeout)
{
    mpiInitializedCheck();
    bool blockingCall = false;
#ifdef USE_BLOCK_SEND
        blockingCall = true;
#endif
    CTimeMon tm(timeout);
    unsigned remaining;
    int target = getRank(dstRank); int tag = getMPITag(mptag);
    bool timedout = false; bool error = false; bool canceled = false; bool completed;
    bool bufferedSendComplete = false;

#ifndef DUAL_MESSAGE_SCHEME
    //TODO: find a better way to send the reply tag
    pushMPTag(mbuf, mbuf.getReplyTag());
#endif
    MPI_Request req;
    while(!bufferedSendComplete && !canceled && !(timedout = tm.timedout(&remaining)))
    {
        if (blockingCall && (target == rank()))
            blockingCall = false;               //incase sending to its self and recv is called by the same thread
        int errorCode;
        ControlMessage cm = ControlMessage(getNextDataTag(),mbuf.getReplyTag());
#ifdef DUAL_MESSAGE_SCHEME
        int dataTag = cm.dataTag;
#else
        int dataTag = tag;
#endif
        if (blockingCall)
        {
#ifdef DUAL_MESSAGE_SCHEME
            errorCode = MPI_Send(&cm, sizeof(ControlMessage), MPI_BYTE, target, tag, this->get());
#endif
            errorCode = MPI_Send(mbuf.bufferBase(), mbuf.length(), MPI_BYTE, target, dataTag, this->get());
        }
        else
        {
#ifdef DUAL_MESSAGE_SCHEME
            errorCode = MPI_Ibsend(&cm, sizeof(ControlMessage), MPI_BYTE, target, tag, this->get(), &req);
#endif
            errorCode = MPI_Ibsend(mbuf.bufferBase(), mbuf.length(), MPI_BYTE, target, dataTag, this->get(), &req);
        }

        int errorClass;
        MPI_Error_class(errorCode, &errorClass);
        if (errorClass == MPI_SUCCESS)
            bufferedSendComplete = true;
        else
        {
            if (errorClass == MPI_ERR_BUFFER)
                usleep(WAIT_DELAY); //retry after giving some time for the buffers to clear up
            else
            {
                throw errorClass;
            }
        }
    }
#ifndef DUAL_MESSAGE_SCHEME
    popMPTag(mbuf);
#endif
    hpcc_mpi::CommStatus status =
            error ? hpcc_mpi::CommStatus::ERROR
                    : (timedout?    hpcc_mpi::CommStatus::TIMEDOUT
                                    :canceled?  hpcc_mpi::CommStatus::CANCELED
                                            :hpcc_mpi::CommStatus::SUCCESS);
    return status;
}

hpcc_mpi::CommStatus hpcc_mpi::MPIComm::readData(rank_t &sourceRank, mptag_t &mptag, CMessageBuffer &mbuf, unsigned timeout)
{
    mpiInitializedCheck();
    bool blockingCall = false;
#ifdef USE_BLOCK_RECV
    blockingCall = true;
#endif
    CTimeMon tm(timeout);
    unsigned remaining;
    bool incomingMessage = false; bool error = false; bool completed = false; bool canceled = false; bool timedout = false;

    tm.timedout(&remaining);
    int source = getRank(sourceRank);
    int tag = getMPITag(mptag);
    CommData* commData = new CommData(source, tag, *this);
    addCommData(commData); //So that it can be cancelled from outside


    MPI_Status stat;
    int errorCode;
    int returnMsgTag;
#ifdef DUAL_MESSAGE_SCHEME
//    if (blockingCall)
//    {
//        //TODO:Allow cancellation
//        errorCode = MPI_Recv(&(commData->cm), sizeof(ControlMessage), MPI_BYTE, source, tag, this->get(), &stat);
//    }
//    else
//    {
    //We are not supporting a blocking receive for the control message since we cannot cancel it at the moment.
    bool notCanceled = commData->lockFromCancellation();
    if (notCanceled)
    {
        errorCode = MPI_Irecv(&(commData->cm), sizeof(ControlMessage), MPI_BYTE, source, tag, this->get(), commData->request());
        commData->releaseCancellationLock();
        stat = waitToComplete(completed, error, canceled, timedout, remaining, commData);
        if (timedout)
            ::cancelComm(commData);
    }
    else
        canceled = true;
//    }
    if (!canceled && !timedout && !error && completed)
    {
        completed = false;
        int dataTag = commData->cm.dataTag;
        returnMsgTag = stat.MPI_TAG;
#else
        int dataTag = tag;
#endif
        stat = hasIncomingData(source, dataTag, *this, incomingMessage, error, canceled, remaining, commData);

        timedout = tm.timedout(&remaining);
        if (incomingMessage && !canceled)
        {
            int size;
            MPI_Get_count(&stat,MPI_BYTE,&size);
            assertex(size>=0);
            mbuf.setLength(size);
            bool notCanceled = commData->lockFromCancellation();
            if (notCanceled)
            {
                if (blockingCall)
                {
                    //TODO:Allow cancellation
                    int errorCode = MPI_Recv(mbuf.bufferBase(), mbuf.length(), MPI_BYTE, source, dataTag, this->get(), &stat);
                    int errorClass;
                    MPI_Error_class(errorCode, &errorClass);
                    completed = (errorClass == MPI_SUCCESS);
                    error = !completed;
                }
                else
                    MPI_Irecv(mbuf.bufferBase(), mbuf.length(), MPI_BYTE, source, dataTag, this->get(), commData->request());

                commData->releaseCancellationLock();
            }
            canceled = !notCanceled;
            tm.timedout(&remaining);

            if (!blockingCall)
                stat = waitToComplete(completed, error, canceled, timedout, remaining, commData);

            if (!canceled)
            {   //if it was canceled by another thread commData would have cleanedup after itself so nothing to do here.
                if (!error && completed)
                {
                    notCanceled = commData->lockFromCancellation();
                    if (notCanceled)
                    {
#ifndef DUAL_MESSAGE_SCHEME
                        returnMsgTag = stat.MPI_TAG;
#endif
                        sourceRank = stat.MPI_SOURCE;
                        mptag = getMPTag(returnMsgTag);
                        //set the reply tag
#ifdef DUAL_MESSAGE_SCHEME
                        mbuf.setReplyTag(commData->cm.replyTag);
#else
                        mbuf.setReplyTag(popMPTag(mbuf));
#endif
                        commData->setCompleted();
                        commData->releaseCancellationLock();
                    }
                    canceled = !notCanceled;
                } else if (timedout)
                {
                    ::cancelComm(commData);
                    canceled = true;
                }
            }
#ifdef DUAL_MESSAGE_SCHEME
        }
#endif
        if (!canceled)
        {
            popCommData(commData);
            delete commData;
        }
    }
    hpcc_mpi::CommStatus status =
            error ? hpcc_mpi::CommStatus::ERROR
                    : (completed? (canceled? hpcc_mpi::CommStatus::CANCELED
                                             : hpcc_mpi::CommStatus::SUCCESS)
                                  : hpcc_mpi::CommStatus::TIMEDOUT);
    return status;
}

void hpcc_mpi::MPIComm::barrier()
{
    mpiInitializedCheck();
    MPI_Barrier(this->get());
}

void hpcc_mpi::mpiInitializedCheck(bool isInitialized)
{
    int flag;
    MPI_Initialized(&flag);
    assert((flag > 0) == isInitialized);
}

void hpcc_mpi::mpiFinalizedCheck(bool isFinalized)
{
    int flag;
    MPI_Finalized(&flag);
    assert((flag > 0) == isFinalized);
}
