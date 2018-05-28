/*##############################################################################

    HPCC SYSTEMS software Copyright (C) 2012 HPCC Systems®.

    Licensed under the Apache License, Version 2.0 (the "License");
    you may not use this file except in compliance with the License.
    You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

    Unless required by applicable law or agreed to in writing, software
    distributed under the License is distributed on an "AS IS" BASIS,
    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
    See the License for the specific language governing permissions and
    limitations under the License.
############################################################################## */

#define mp_decl DECL_EXPORT
#include "platform.h"
#include "jlib.hpp"
#include "jlog.hpp"

#include "mpbase.hpp"
#include "mpi_wrapper.hpp"

static INode *MyNode=NULL;
static INode *NullNode=NULL;

class MPIGroup: implements IGroup, public CInterface{
protected: friend class CNodeIterator;
    rank_t count;
    mutable rank_t myrank;
public:
    IMPLEMENT_IINTERFACE;

    MPIGroup(){
        count = hpcc_mpi::size();
        myrank = hpcc_mpi::rank();
    }

    ~MPIGroup(){
    }

    rank_t ordinality()  const{ 
        return count; 
    }
    
    rank_t rank(const SocketEndpoint &ep) const {
        UNIMPLEMENTED;
    }
    
    rank_t rank(INode *node) const  { 
        UNIMPLEMENTED;
    }
    
    rank_t rank() const { 
        return myrank;
    }

    GroupRelation compare(const IGroup *grp) const {
        UNIMPLEMENTED;
    }

    bool equals(const IGroup *grp) const{
        UNIMPLEMENTED;
    }

    void translate(const IGroup *othergroup, rank_t nranks, const rank_t *otherranks, rank_t *resranks ) const {
        UNIMPLEMENTED;
    }

    IGroup *subset(rank_t start,rank_t num) const{
        UNIMPLEMENTED;
    }

    virtual IGroup *subset(const rank_t *order,rank_t num) const{
        UNIMPLEMENTED;
    }

    virtual IGroup *combine(const IGroup *grp) const{
        UNIMPLEMENTED;
    }

    bool isMember(INode *node) const{
        UNIMPLEMENTED;
    }

    bool isMember() const{
        UNIMPLEMENTED;
    }

    unsigned distance(const IpAddress &ip) const{
        UNIMPLEMENTED;
    }

    unsigned distance(const IGroup *grp) const{
        UNIMPLEMENTED;
    }

    IGroup *diff(const IGroup *g) const{
        UNIMPLEMENTED;
    }

    bool overlaps(const IGroup *grp) const{
        UNIMPLEMENTED;        
    }


    IGroup *intersect(const IGroup *g) const{
        UNIMPLEMENTED;
    }

    IGroup *swap(rank_t r1,rank_t r2) const{
        UNIMPLEMENTED;
    }

    virtual IGroup *add(INode *node) const {
        UNIMPLEMENTED;
    }

    virtual IGroup *add(INode *node,unsigned pos) const {
        UNIMPLEMENTED;
    }

    virtual IGroup *remove(unsigned pos) const{
        UNIMPLEMENTED;
    }

    virtual IGroup *rotate(int num) const{
        UNIMPLEMENTED;
    }


    INode &queryNode(rank_t r) const {
        UNIMPLEMENTED;
    }
    INode *getNode(rank_t r) const {
        UNIMPLEMENTED;
    }

    StringBuffer &getText(StringBuffer &text) const {
        UNIMPLEMENTED;
    }

    static IGroup *fromText(const char *s,unsigned defport) {
        UNIMPLEMENTED;
    }

    void serialize(MemoryBuffer &tgt) const{
        UNIMPLEMENTED;
    }
    
    static IGroup *deserialize(MemoryBuffer &src) {
        UNIMPLEMENTED;
    }

    void getSocketEndpoints(SocketEndpointArray &sea) const{
        UNIMPLEMENTED;
    }
    INodeIterator *getIterator(rank_t start=0,rank_t num=RANK_NULL) const {
        UNIMPLEMENTED; 
    }
    
};

IGroup *createIGroup(rank_t num,INode **nodes){
    return new MPIGroup();
}

IGroup *createIGroup(rank_t num,const SocketEndpoint *ep){
    UNIMPLEMENTED; 
}

IGroup *createIGroup(SocketEndpointArray &epa){
    UNIMPLEMENTED; 
}

IGroup *createIGroup(const char *endpointlist,unsigned short defport){
    UNIMPLEMENTED; 
}

class MPNode: implements INode, public CInterface
{
protected: friend class MPNodeCache;
    SocketEndpoint ep;
public:
    IMPLEMENT_IINTERFACE;
    MPNode(const SocketEndpoint &_ep)
        : ep(_ep)
    {
#ifdef _DEBUG
//      assertex(!_ep.LoopBack()); // localhost not supported for nodes
#endif
    }
    bool equals(const INode *n) const { return endpoint().equals(n->endpoint()); }
    void serialize(MemoryBuffer &tgt)
    {
        ep.serialize(tgt);
    }
    static MPNode *deserialize(MemoryBuffer &src);
    unsigned getHash() const { return ep.hash(0); }
    virtual bool isHost() const
    {
        return ep.isHost();
    }
    virtual bool isLocalTo(INode *to) const
    {
        return ep.ipequals(to->endpoint());
    }
    const SocketEndpoint &endpoint() const { return ep; }

};

class MPNodeCache: public SuperHashTableOf<MPNode,SocketEndpoint>
{
    CriticalSection sect;
public:
    ~MPNodeCache()
    {
        _releaseAll();
    }

    void onAdd(void *)
    {
        // not used
    }

    void onRemove(void *e)
    {
        MPNode &elem=*(MPNode *)e;      
        elem.Release();
    }

    unsigned getHashFromElement(const void *e) const
    {
        const MPNode &elem=*(const MPNode *)e;      
        return elem.ep.hash(0);
    }

    unsigned getHashFromFindParam(const void *fp) const
    {
        return ((const SocketEndpoint *)fp)->hash(0);
    }

    const void * getFindParam(const void *p) const
    {
        const MPNode &elem=*(const MPNode *)p;      
        return &elem.ep;
    }

    bool matchesFindParam(const void * et, const void *fp, unsigned) const
    {
        return ((MPNode *)et)->ep.equals(*(SocketEndpoint *)fp);
    }

    IMPLEMENT_SUPERHASHTABLEOF_REF_FIND(MPNode,SocketEndpoint);

    MPNode *lookup(const SocketEndpoint &ep)
    {
        CriticalBlock block(sect);
        MPNode *item=SuperHashTableOf<MPNode,SocketEndpoint>::find(&ep);
        if (!item) {
            item = new MPNode(ep);
            add(*item);
        }
        return LINK(item);
    }


} *NodeCache = NULL;

MPNode *MPNode::deserialize(MemoryBuffer &src)
{
    SocketEndpoint ep;
    ep.deserialize(src);
    if (NodeCache)
        return NodeCache->lookup(ep);
    return new MPNode(ep);
}


class CNodeIterator : implements INodeIterator, public CInterface
{
    Linked<IGroup> parent;
    rank_t start;
    rank_t num;
    rank_t pos;
public:
    IMPLEMENT_IINTERFACE;
    CNodeIterator(IGroup *_parent,rank_t _start,rank_t _num)
        :   parent(_parent)
    {
        start = _start;
        num = _num;
        if ((num==RANK_NULL)||(num+start>parent->ordinality()))
            num = (start<parent->ordinality())?parent->ordinality()-start:0;
    }
    bool first()
    {
        pos = 0;
        return (pos<num);
    }
    bool next()
    {
        pos++;
        return (pos<num);
    }
    bool isValid()
    {
        return (pos<num);
    }
    INode &query()
    {
        return parent->queryNode(start+pos);
    }
    INode &get()
    {
        return *parent->getNode(start+pos);
    }
};

INodeIterator *IGroup::getIterator(rank_t start,rank_t num) const{
    UNIMPLEMENTED; 
}



INode *deserializeINode(MemoryBuffer &src)
{
    return MPNode::deserialize(src);
}

INode *createINode(const SocketEndpoint &ep)
{
    if (NodeCache)
        return NodeCache->lookup(ep);
    return new MPNode(ep);
}

INode *createINodeIP(const IpAddress &ip,unsigned short port=0)
{
    SocketEndpoint ep(port,ip);
    return createINode(ep);
}

INode *createINode(const char *name,unsigned short port)
{
    SocketEndpoint ep(name,port);
    return createINode(ep);
}

IGroup *deserializeIGroup(MemoryBuffer &src)
{
    return MPIGroup::deserialize(src);
}



void initMyNode(unsigned short port)
{
    setNodeCaching(port != 0);
    ::Release(MyNode);
    MyNode = NULL;
    if (port) {
        SocketEndpoint ep(port);
        MyNode = new MPNode(ep);
        if (ep.isLoopBack()) {
            DBGLOG("MP Warning - localhost used for MP host address, NIC adaptor not identified");
        }
        queryNullNode();
    }
    else
    {
        ::Release(NullNode);
        NullNode = NULL;
    }
}

INode *queryMyNode()
{
    return MyNode;
}

INode *queryNullNode()
{
    if (!NullNode) {
        SocketEndpoint ep;
        NullNode = new MPNode(ep);
    }
    return NullNode;
}

void setNodeCaching(bool on)
{
    if (on) {
        if (!NodeCache)
            NodeCache = new MPNodeCache();
    }
    else { 
        MPNodeCache *nc = NodeCache;
        NodeCache = NULL;
        delete nc;
    }
}

