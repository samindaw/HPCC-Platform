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

//----------Functions and Data structures managing communication data related to Send/Recv Communications in orogress-----------//

// Data structure to keep the data relating to send/receive communications
class CommData
{
private:
    bool cancellationLock = false;
    bool cancellationInProgress = false;
    MPI_Request req;           // persistent request object to keep track of ongoing MPI call
    const hpcc_mpi::MPIComm& comm;                // MPI communicator

    CriticalSection dataChangeLock;
    void lock(){dataChangeLock.enter();}
    void unlock(){dataChangeLock.leave();}

public:    
    int rank;                       // source/destination rank of the processor
    int tag;                        // MPI tag information

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
    void notifyCancellation()
    {
        lock();
        _TF("(rank, tag)",rank, tag);
        cancellationInProgress = true;
        unlock();
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
    _TF("_popCommData", index);
    CommData* ret = NULL;
    assertex(index>=0);
    assertex(index < asyncCommData.size());
    ret = asyncCommData[index];
    asyncCommData.erase(asyncCommData.begin() + index);
    return ret;
}

void addCommData(CommData *commData)
{
    _TF("addCommData", commData->rank, commData->tag);
    CriticalBlock block(commDataLock);

    //TODO Do a cleanup while we are at it
//    int size = asyncCommData.size(); int completed; MPI_Status stat;
//    for(int i=(size-1); i>=0 ; i--)
//    {
//        if (!(asyncCommData[i]->probingProgress))
//        {
//            completed = 0;
//            assertex(asyncCommData[i]->request != NULL);
//            _T("asyncCommData[i]->request="<<*(asyncCommData[i]->request)<<" mem_address="<<asyncCommData[i]->request);
//            _T("send="<<asyncCommData[i]->isSend()<<" rank="<<asyncCommData[i]->rank<<" tag="<<asyncCommData[i]->tag);
//            bool error = (MPI_Test(asyncCommData[i]->request, &completed, &stat)!= MPI_SUCCESS);
//            if (completed || error) //unlikely an error would occur
//            {
//                delete _popCommData(i);
//            }
//        }
//    }

    asyncCommData.push_back(commData);
}

void popCommData(int rank, int tag, const hpcc_mpi::MPIComm& comm, std::vector<CommData*> &matchedResults)
{
   _TF("popCommData", rank, tag);
   CriticalBlock block(commDataLock);
   size_t size = asyncCommData.size();
   for(size_t i=(size-1); i>=0; i--)
   {
       if (asyncCommData[i]->isEqual(rank, tag))
           matchedResults.push_back(_popCommData(i));
   }
}

CommData* popCommData(CommData * commData)
{
    _TF("popCommData(CommData * commData)", commData->rank, commData->tag);
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

#define TAG_UB 1000000 //INT_MAX //MPI_TAG_UB
int getRank(rank_t sourceRank)
{
    _TF("getRank", sourceRank);
    if (sourceRank == RANK_ALL)
        return MPI_ANY_SOURCE;
    else
        return sourceRank;
}

#define SPECIAL_TAG_BASE (TAG_UB/2) // reserve half range for special tags (inc. reply tags)
int getMPITag(mptag_t mptag)
{
    _TF("getMPITag", mptag);

    if (mptag == TAG_ALL)
        return MPI_ANY_TAG;
    else
    {
        int tag = (int)mptag;
        if (tag < 0)
        {
            assertex(-tag < TAG_UB-SPECIAL_TAG_BASE); // check not out of range
            return -tag + SPECIAL_TAG_BASE;
        }
        else
        {
            assertex(tag <= TAG_UB);
            return tag;
        }
    }
}

mptag_t getMPTag(int tag)
{
    _TF("getMPTag", tag);
    unsigned returnVal;
    if (tag == MPI_ANY_TAG)
        returnVal = TAG_ALL;
    else
    {
        if (tag > SPECIAL_TAG_BASE)
            returnVal = (mptag_t) -(tag-SPECIAL_TAG_BASE);
        else
            returnVal = tag;
    }
    return (mptag_t)returnVal;
}

MPI_Status waitToComplete( bool& completed, bool& error, bool& canceled, bool& timedout, unsigned timeout, CommData *commData)
{
    _TF("waitToComplete", completed, error, canceled, timeout);
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
    _TF("hasIncomingData", sourceRank, mptag, incomingMessage, error, timeout);

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
    _TF("hasIncomingMessage", sourceRank, mptag);
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
    _TF("cancelComm", send, rank, mptag);
    mpiInitializedCheck();
    int r = getRank(rank);
    int tag = getMPITag(mptag);
    std::vector<CommData*> commDataList;
    popCommData(r, tag, *this, commDataList);
    bool success = true;
    for (auto &commData: commDataList)
    {
        commData->notifyCancellation();              // In case sendData and readData functions still in progress we want to let
                                                                    // them know that we are about to screw up their plans
        if (!::cancelComm(commData))
            success = false;                        // If we managed to cancel everything then the cancellation was successful
        delete commData;
    }
    return success;
}

rank_t hpcc_mpi::MPIComm::rank()
{
    _TF("rank");
    mpiInitializedCheck();
    int rank;
    MPI_Comm_rank(this->get(), &rank);
    return rank;
}

rank_t hpcc_mpi::MPIComm::size()
{
    _TF("size");
    mpiInitializedCheck();
    int size;
    MPI_Comm_size(this->get(), &size);
    return size;
}

#define MPI_BUFFER_SIZE 16777216        // 16MB
void* mpi_buffer_data;

void hpcc_mpi::initialize()
{
    _TF("initialize");
    mpiInitializedCheck(false);

    int required = MPI_THREAD_MULTIPLE;
    int provided;

    MPI_Init_thread(NULL, NULL, required, &provided);
    assertex(provided == required);

    //Allocate memory for the MPI buffered send
    mpi_buffer_data = malloc(MPI_BUFFER_SIZE);
    MPI_Buffer_attach(mpi_buffer_data, MPI_BUFFER_SIZE);

#ifdef DEBUG
    MPI_Comm_rank(MPI_COMM_WORLD, &global_proc_rank);
#endif
}

void hpcc_mpi::MPIComm::setErrorHandler(MPI_Errhandler handler)
{
    mpiInitializedCheck();
    MPI_Comm_set_errhandler(this->get(),handler);
}

void hpcc_mpi::finalize()
{
    _TF("finalize");
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
    _TF("sendData", dstRank, mptag, mbuf.getReplyTag(), timeout);
    mpiInitializedCheck();
    CTimeMon tm(timeout);
    unsigned remaining;
    int target = getRank(dstRank); int tag = getMPITag(mptag);
    bool timedout = false; bool error = false; bool canceled = false; bool completed;
    bool bufferedSendComplete = false;

    //TODO: find a better way to send the reply tag
    pushMPTag(mbuf, mbuf.getReplyTag());
    MPI_Request req;
    while(!bufferedSendComplete && !canceled && !(timedout = tm.timedout(&remaining)))
    {
        int error = MPI_Ibsend(mbuf.bufferBase(), mbuf.length(), MPI_BYTE, target, tag, this->get(), &req);
        int errorClass;
        MPI_Error_class(error, &errorClass);
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
    popMPTag(mbuf);

    hpcc_mpi::CommStatus status =
            error ? hpcc_mpi::CommStatus::ERROR
                    : (timedout?    hpcc_mpi::CommStatus::TIMEDOUT
                                    :canceled?  hpcc_mpi::CommStatus::CANCELED
                                            :hpcc_mpi::CommStatus::SUCCESS);
    return status;
}

hpcc_mpi::CommStatus hpcc_mpi::MPIComm::readData(rank_t &sourceRank, mptag_t &mptag, CMessageBuffer &mbuf, unsigned timeout)
{
    _TF("readData", sourceRank, mptag, timeout);
    mpiInitializedCheck();
    CTimeMon tm(timeout);
    unsigned remaining;
    bool incomingMessage = false; bool error = false; bool completed = false; bool canceled = false; bool timedout = false;

    tm.timedout(&remaining);
    int source = getRank(sourceRank);
    int tag = getMPITag(mptag);

    CommData* commData = new CommData(source, tag, *this);
    addCommData(commData); //So that it can be cancelled from outside

    MPI_Status stat = hasIncomingData(source, tag, *this, incomingMessage, error, canceled, remaining, commData);

    timedout = tm.timedout(&remaining);
    if (incomingMessage && !canceled)
    {
        int size;
        MPI_Get_count(&stat,MPI_BYTE,&size);
        assertex(size>0);
        mbuf.setLength(size);
        bool notCanceled = commData->lockFromCancellation();
        if (notCanceled)
        {
            MPI_Irecv(mbuf.bufferBase(), mbuf.length(), MPI_BYTE, source, tag, this->get(), commData->request());
            commData->releaseCancellationLock();
        }
        canceled = !notCanceled;
        tm.timedout(&remaining);

        MPI_Status stat = waitToComplete(completed, error, canceled, timedout, remaining, commData);
        if (!canceled)
        {   //if it was canceled by another thread commData would have cleanedup after itself so nothing to do here.
            if (!error && completed)
            {
                notCanceled = commData->lockFromCancellation();
                if (notCanceled)
                {
                    sourceRank = stat.MPI_SOURCE;
                    mptag = getMPTag(stat.MPI_TAG);
                    //pop the reply tag
                    mbuf.setReplyTag(popMPTag(mbuf));
                    commData->releaseCancellationLock();
                }
                canceled = !notCanceled;
            } else if (timedout)
            {
                ::cancelComm(commData);
                canceled = true;
            }
        }

    }
    if (!canceled)
    {
        popCommData(commData);
        delete commData;
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
    _TF("barrier");
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
