/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include <mpi/mpi.h>
#include "mpi_wrapper.hpp"
#include <cstdlib>
#include <queue>
#include <typeinfo>

#include "mputil.hpp"

#define WAIT_DELAY 100
#define PROBING_DELAY 100

//----------Functions and Data structures managing communication data related to Send/Recv Communications in orogress-----------//

// Data structure to keep the data relating to send/receive communications
class CommData
{
private:
    bool send;                      // TRUE => relates to a send communication | FALSE => relates to receive communication
    bool locked = false;
    bool cancellationLock = false;
    bool cancellationInProgress = false;
    CMessageBuffer* data = NULL;              // Data structure which points to the sent/recv buffer
    CriticalSection dataChangeLock;
    void lock(){if (!locked) dataChangeLock.enter(); locked=true;}
    void unlock(){if (locked) dataChangeLock.leave(); locked=false;}
public:    
    bool isSend(){return send;}
    bool isReceive(){return !isSend();}
    bool isEqual(bool _send, int _rank, int _tag, MPI_Comm _comm) {return (send==_send) && (_rank==MPI_ANY_SOURCE || rank==_rank) && (_tag==MPI_ANY_TAG || tag==_tag) && (comm == _comm);}

    int rank;                       // source/destination rank of the processor
    int tag;                        // MPI tag information
    MPI_Request *request;           // persistent request object to keep track of ongoing MPI call
    MPI_Comm comm;                  // MPI communicator


    CommData(bool _send, int _rank, int _tag, MPI_Comm _comm):
        send(_send), rank(_rank), tag(_tag), comm(_comm)
    {
        request = new MPI_Request();
    }

    CommData(int _rank, int _tag, MPI_Comm _comm): rank(_rank), tag(_tag), comm(_comm)
    {
        request = new MPI_Request();
        send = false;
    }

    void loadBuffer(CMessageBuffer &buf)
    {
        data  = &buf;
    }

    void* buffer()
    {
        return data->bufferBase();
    }
    void updateBufferLength(size_t length)
    {
        data->setLength(length);
    }
    size_t bufferLength()
    {
        return data->length();
    }

    void notifyCancellation()
    {
        lock();
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

    ~CommData()
    {
        if (request) delete request;
    }
};

std::vector<CommData*> asyncCommData; // CommData list to manage while send/recv communication in progress
CriticalSection commDataLock;         // A mutex lock for the index list above

CommData* _popCommData(int index)
{
    _TF("_popCommData", index);
    CommData* ret = NULL;
    if (index != -1) {
        ret = asyncCommData[index];
        asyncCommData.erase(asyncCommData.begin() + index);
    }
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

CommData* popCommData(int index)
{
    _TF("popCommData", index);
    CriticalBlock block(commDataLock);
    return _popCommData(index);
}

std::vector<CommData*> popCommData(bool send, int rank, int tag, MPI_Comm comm)
{
   _TF("popCommData", rank, tag);
   CriticalBlock block(commDataLock);
   std::vector<CommData*> matchedResults = std::vector<CommData*>();
   int size = asyncCommData.size();
   for(int i=(size-1); i>=0; i--)
   {
       if (asyncCommData[i]->isEqual(send, rank, tag, comm))
           matchedResults.push_back(_popCommData(i));
   }
   return matchedResults;
}

CommData* popCommData(CommData * commData)
{
    _TF("popCommData(CommData * commData)", commData->rank, commData->tag, commData->comm);
    CriticalBlock block(commDataLock);
    int size = asyncCommData.size();
    for(int i=(size-1); i>=0; i--)
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

#define TAG_UB 100000 //INT_MAX //MPI_TAG_UB
int getRank(rank_t sourceRank)
{
    _TF("getRank", sourceRank);
    if (sourceRank == RANK_ALL)
        return MPI_ANY_SOURCE;
    else
        return sourceRank;
}

int getTag(mptag_t mptag)
{
    _TF("getTag", mptag);
    int returnVal;
    unsigned tag = mptag;
    if (mptag == TAG_ALL)
        returnVal = MPI_ANY_TAG;
    else
    {
        if (tag > TAG_UB)
        {
            tag= TAG_UB-((unsigned)-tag);                         //MPI doesn't allow custom tags with negative values or beyond MPI_TAG_UB. Thus we shift it to the MPI_TAG_UB range
        }
        returnVal = static_cast<int>(tag);

    }
    return returnVal;
}

mptag_t getTag(int tag)
{
    _TF("getTag", tag);
    unsigned returnVal;
    if (tag == MPI_ANY_TAG)
        returnVal = TAG_ALL;
    else
    {
        if (tag > (TAG_UB/2))
        {
            returnVal = (unsigned) - (TAG_UB-tag);                         //MPI doesn't allow custom tags with negative values
        } else
            returnVal = tag;

    }
    return (mptag_t)returnVal;
}

MPI_Status waitToComplete( bool& completed, bool& error, bool& canceled, bool& timedout, unsigned timeout, CommData *commData)
{
    _TF("mpi_wrapper:waitToComplete", completed, error, canceled, timeout);
    CTimeMon tm(timeout);
    MPI_Status stat;
    int flag;
    unsigned remaining;
    bool noCancellation;
    while ((noCancellation = commData->lockFromCancellation()) && !(completed || error || (timedout = tm.timedout(&remaining))))
    {
        error = (MPI_Test(commData->request, &flag, &stat) != MPI_SUCCESS);
        completed = (flag > 0);
        commData->releaseCancellationLock();
        usleep(WAIT_DELAY);
    }
    if (completed)
    {
        MPI_Test_cancelled(&stat, &flag);
        canceled = (flag > 0);
    }else
    {
        canceled = !noCancellation;
    }

    return stat;
}

MPI_Status hasIncomingData(int sourceRank, int mptag, MPI_Comm comm,
        bool& incomingMessage, bool& error, bool& canceled, unsigned timeout, CommData *commData)
{
    _TF("mpi_wrapper:hasIncomingData", sourceRank, mptag, incomingMessage, error, timeout);

    CTimeMon tm(timeout);
    MPI_Status stat;

    int flag;
    unsigned remaining;
    bool noCancellation;

    while ((noCancellation = commData->lockFromCancellation()) && !(incomingMessage || error || (timeout !=0 && tm.timedout(&remaining))))
    {
        error = (MPI_Iprobe(sourceRank, mptag, comm, &flag, &stat) != MPI_SUCCESS);
        incomingMessage = (flag > 0);
        commData->releaseCancellationLock();
        if (timeout == 0) break;
        usleep(PROBING_DELAY);
    }
    canceled = !noCancellation;
    return stat;
}

//----------------------------------------------------------------------------//

/** See mpi_wrapper.hpp header file for function descriptions of the following **/

bool hpcc_mpi::hasIncomingMessage(rank_t &sourceRank, mptag_t &mptag, MPI_Comm comm)
{
    _TF("mpi_wrapper:hasIncomingMessage", sourceRank, mptag);
    bool incomingMessage = false; bool error = false; bool canceled = false;
    int source = getRank(sourceRank);
    int tag = getTag(mptag);
    CommData* tmpCommData = new CommData(source, tag, comm);
    MPI_Status stat = hasIncomingData(source, tag, comm, incomingMessage, error, canceled, 0, tmpCommData);
    if (incomingMessage)
    {
        sourceRank = stat.MPI_SOURCE;
        mptag = mptag_t(stat.MPI_TAG);
    }

    return incomingMessage;
}

int cancelComm(CommData* commData)
{
    int ret = true;
    if (commData)
    {
        MPI_Status stat;
        int completed;
        bool error = (MPI_Test(commData->request, &completed, &stat)
                != MPI_SUCCESS);
        if (!error && !completed) {
            error = (MPI_Cancel(commData->request) != MPI_SUCCESS);
            if (!error) {
                // MPI framework managed to successfully cancel
                MPI_Request_free(commData->request);
            }
        }
        ret = !error;
    }

    return ret;
}

bool hpcc_mpi::cancelComm(bool send, rank_t rank, mptag_t mptag, MPI_Comm comm)
{
    _TF("mpi_wrapper:cancelComm", send, rank, mptag, comm);
    int r = getRank(rank);
    int tag = getTag(mptag);
    std::vector<CommData*> commDataList = popCommData(send, r, tag, comm);
    bool success = true;
    for(int i=0; i<commDataList.size(); i++)
    {
        commDataList[i]->notifyCancellation();              // In case sendData and readData functions still in progress we want to let
                                                            // them know that we are about to screw up their plans
        success = success && cancelComm(commDataList[i]);   // If we managed to cancel everything then the cancellation was successful
        usleep(WAIT_DELAY);                                 // Wait for a short while for active threads to exit send/recv function
        delete commDataList[i];
    }
    return success;
}

rank_t hpcc_mpi::rank(MPI_Comm comm)
{
    _TF("mpi_wrapper:rank");
    int rank;
    MPI_Comm_rank(comm, &rank);
    return rank;
}

rank_t hpcc_mpi::size(MPI_Comm comm)
{
    _TF("mpi_wrapper:size");
    int size;
    MPI_Comm_size(comm, &size);
    return size;
}

void hpcc_mpi::initialize(bool withMultithreading)
{
    _TF("mpi_wrapper:initialize", withMultithreading);
    if (withMultithreading)
    {
        int required = MPI_THREAD_MULTIPLE;
        int provided;
        MPI_Init_thread(NULL,NULL, required, &provided);
        assertex(provided == required);
    }else
    {
        MPI_Init(NULL,NULL);
    }
#ifdef DEBUG
    MPI_Comm_rank(MPI_COMM_WORLD, &global_proc_rank);
#endif
}

void hpcc_mpi::finalize()
{
    _TF("mpi_wrapper:finalize");
    MPI_Finalize();
}

hpcc_mpi::CommStatus hpcc_mpi::sendData(rank_t dstRank, mptag_t mptag, CMessageBuffer &mbuf, MPI_Comm comm, unsigned timeout)
{
    _TF("mpi_wrapper:sendData", dstRank, mptag, timeout);
    CTimeMon tm(timeout);
    unsigned remaining;
    int target = getRank(dstRank); int tag = getTag(mptag);

    CommData* commData = new CommData(true, target, tag, comm);
    commData->loadBuffer(mbuf);
    bool timedout = false; bool error = false; bool canceled = false; bool completed;

    error  = (MPI_Isend(commData->buffer(), commData->bufferLength(), MPI_BYTE, target, tag, comm, commData->request) != MPI_SUCCESS);
    timedout = tm.timedout(&remaining);

    addCommData(commData); //So that it can be cancelled from outside

    if (timeout == MP_WAIT_FOREVER) //if blocking send is requested
    {
        waitToComplete(completed, error, canceled, timedout, timeout, commData);
        if (completed)
        {
            popCommData(commData);
            delete commData;
        }
    }
    if (!error && timedout)
    {
        popCommData(commData);
        cancelComm(commData);
    }

    hpcc_mpi::CommStatus status =
            error ? hpcc_mpi::CommStatus::ERROR
                    : (timedout?    hpcc_mpi::CommStatus::TIMEDOUT
                                    :canceled?  hpcc_mpi::CommStatus::CANCELED
                                            :hpcc_mpi::CommStatus::SUCCESS);
    return status;
}

hpcc_mpi::CommStatus hpcc_mpi::readData(rank_t sourceRank, mptag_t mptag, CMessageBuffer &mbuf, MPI_Comm comm, unsigned timeout)
{
    _TF("mpi_wrapper:readData", sourceRank, mptag, timeout);
    CTimeMon tm(timeout);
    unsigned remaining;
    bool incomingMessage = false; bool error = false; bool completed = false; bool canceled = false; bool timedout = false;

    tm.timedout(&remaining);
    int source = getRank(sourceRank);
    int tag = getTag(mptag);

    CommData* commData = new CommData(source, tag, comm);
    addCommData(commData); //So that it can be cancelled from outside

    MPI_Status stat = hasIncomingData(source, tag, comm, incomingMessage, error, canceled, remaining, commData);

    timedout = tm.timedout(&remaining);
    if (incomingMessage && !canceled)
    {
        int size;
        MPI_Get_count(&stat, MPI_BYTE, &size);
        assertex(size>0);
        commData->loadBuffer(mbuf);
        commData->updateBufferLength(size);
        error  = (MPI_Irecv(commData->buffer(), commData->bufferLength(), MPI_BYTE, source, tag, comm, commData->request) != MPI_SUCCESS);
        tm.timedout(&remaining);

        //-----------------------------------------Blocking Recv----------------------------------------//
        MPI_Status stat = waitToComplete(completed, error, canceled, timedout, remaining, commData);
        if (!canceled)
        {   //if it was canceled by another thread commData would have cleanedup after itself so nothing to do here.
            if (!error && completed)
            {
                bool noCancellation = commData->lockFromCancellation();
                if (noCancellation)
                {
                    SocketEndpoint ep(stat.MPI_SOURCE);
                    mbuf.init(ep, (mptag_t)getTag(stat.MPI_TAG), TAG_REPLY_BASE);
                    commData->releaseCancellationLock();
                }
                canceled = !noCancellation;
            } else if (timedout)
            {
                cancelComm(commData);
                canceled = true;
            }
        }
        //--------------------------------------END: Blocking Recv--------------------------------------//

        /*-----------------We will not support a scenario for a asynchronous recv-------------------
         *
                if (timeout == MP_ASYNC_SEND)
                {
                    commData->probingProgress = false;
                }else
                {
                    //-----------------------------------------Blocking Recv----------------------------------------//
                    MPI_Status stat = waitToComplete(commData->request, completed, error, canceled, timedout, remaining);
                    if (!canceled)
                    {   //if it was canceled by another thread commData would have cleanedup after itself so nothing to do here.
                        _T("Irecv completed="<<completed<<" error="<<error<<" canceled="<<canceled);
                        popCommData(commData);
                        if (!error && completed)
                        {
                            mbuf.reset();
                            mbuf.append(size,commData->data);
                            SocketEndpoint ep(stat.MPI_SOURCE);
                            mbuf.init(ep, (mptag_t)(stat.MPI_TAG), TAG_REPLY_BASE);
                        } else if (timedout)
                        {
                            cancelComm(commData);
                        }
                        delete commData;
                    }
                    //--------------------------------------END: Blocking Recv--------------------------------------//
                }

         */
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

void hpcc_mpi::barrier(MPI_Comm comm)
{
    _TF("mpi_wrapper:barrier");
    MPI_Barrier(comm);
}