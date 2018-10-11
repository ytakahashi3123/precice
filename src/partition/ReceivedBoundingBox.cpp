#include "partition/ReceivedBoundingBox.hpp"
#include "com/CommunicateBoundingBox.hpp"
#include "com/CommunicateMesh.hpp"
#include "com/Communication.hpp"
#include "utils/MasterSlave.hpp"
#include "m2n/M2N.hpp"
#include "utils/EventTimings.hpp"
#include "mapping/Mapping.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/Vertex.hpp"
#include "mesh/Edge.hpp"
#include "mesh/Triangle.hpp"
#include "utils/Helpers.hpp"
#include "utils/MasterSlave.hpp"
#include <vector>
#include <map>

using precice::utils::Event;

namespace precice {
namespace partition {

logging::Logger ReceivedBoundingBox:: _log ( "precice::partition::ReceivedBoundingBox" );

ReceivedBoundingBox::ReceivedBoundingBox
(
  mesh::PtrMesh mesh, double safetyFactor, GeometricFilter geometricFilter)
:
  Partition (mesh),
  //_bb(mesh->getBoundingBox()),
  _bb(mesh->getDimensions(), std::make_pair(std::numeric_limits<double>::max(),
                                            std::numeric_limits<double>::lowest())),
  _dimensions(mesh->getDimensions()),
  _safetyFactor(safetyFactor),
  _geometricFilter(geometricFilter)
{}

void ReceivedBoundingBox::communicateBoundingBox()
{
  TRACE();
  Event e("receive global bounding box");
  prepareBoundingBox();
  if (not utils::MasterSlave::_slaveMode) {
    remoteParComSize=0;
    _m2n->getMasterCommunication()->receive(remoteParComSize, 0);
    for (int remoteRank = 0; remoteRank < remoteParComSize; remoteRank++ ) {
      _globalBB[remoteRank]= _bb;
    }
// master receives global bb from other master    
    com::CommunicateBoundingBox(_m2n->getMasterCommunication()).receiveBoundingBoxMap(_globalBB, 0 );
  }
}

void ReceivedBoundingBox::computeBoundingBox()
{
  // handle coupling mode first (i.e. serial participant)
  if (not utils::MasterSlave::_slaveMode && not utils::MasterSlave::_masterMode)
  { //coupling mode
  //some code here
  }

  Event e1("broadcast received bounding box");
   

    if (utils::MasterSlave::_slaveMode)
    {
      utils::MasterSlave::_communication->broadcast(remoteParComSize, 0);

      // initializing the _globalBB
      for (int remoteRank = 0; remoteRank < remoteParComSize; remoteRank++ ) {
      _globalBB[remoteRank]= _bb;
      }
      // receive _globalBB from master
      com::CommunicateBoundingBox(utils::MasterSlave::_communication).broadcastReceiveBoundingBoxMap(_globalBB);
      numberOfRemoteRanks = _globalBB.size();      

      // received bounding boxes are comapred with each rank bb
      if (numberOfRemoteRanks>0) {
        for (auto &other_rank: _globalBB)
        {
          if (CompareBoundingBox(_bb,other_rank.second)) {
            feedback.push_back(other_rank.first);            
          }
        }
      } 
     //send feedback to master
       utils::MasterSlave::_communication->send(feedback, 0);     
    }
    else if (utils::MasterSlave::_masterMode)
    { // Master

      assertion(utils::MasterSlave::_rank==0);
      assertion(utils::MasterSlave::_size>1);
      _m2n->getMasterCommunication()->send(utils::MasterSlave::_size , 0);
      utils::MasterSlave::_communication->broadcast(remoteParComSize);
      com::CommunicateBoundingBox(utils::MasterSlave::_communication).broadcastSendBoundingBoxMap(_globalBB);
      numberOfRemoteRanks = _globalBB.size();      
      
      if (numberOfRemoteRanks >0) {           
        for (auto &other_rank: _globalBB)
        {
          if (CompareBoundingBox(_bb,other_rank.second)) {
            feedback.push_back(other_rank.first);
          }          
        }
        feedbackMap[0]=feedback;
      }
      

      for (int rank_slave=1; rank_slave < utils::MasterSlave::_size ; rank_slave++) {
        utils::MasterSlave::_communication->receive(feedback, rank_slave);
        feedbackMap[rank_slave]=feedback;        
      }
        com::CommunicateBoundingBox(_m2n->getMasterCommunication()).sendFeedbackMap(feedbackMap,0);             
    }
    e1.stop();


    if (utils::MasterSlave::_masterMode)
    {
      for (int i = 0; i < remoteParComSize; i++) {
        int counter=0;
        _m2n->getMasterCommunication()->receive(counter, 0);
        vertexCounters.push_back(counter);
      }   
      utils::MasterSlave::_communication->broadcast(vertexCounters);
    }
    else if (utils::MasterSlave::_slaveMode)
    {    
      utils::MasterSlave::_communication->broadcast(vertexCounters, 0);  
    }  

}

void ReceivedBoundingBox::communicate()
{
  _m2n->broadcastReceiveLocalMesh(*_mesh);
}

void ReceivedBoundingBox::compute()
{
  TRACE(_geometricFilter);

  if (not utils::MasterSlave::_slaveMode) {
    CHECK(_fromMapping.use_count() > 0 || _toMapping.use_count() > 0,
          "The received mesh " << _mesh->getName()
          << " needs a mapping, either from it, to it, or both. Maybe you don't want to receive this mesh at all?")
  }


  // To understand the following steps, it is recommended to look at BU's thesis, especially Figure 69 on page 89 
  // for RBF-based filtering. https://mediatum.ub.tum.de/doc/1320661/document.pdf


  // (0) set global number of vertices before filtering

  // (1) Bounding-Box-Filter
  

  INFO("Broadcast mesh " << _mesh->getName());
    
  if (_geometricFilter == BROADCAST_FILTER) {

    INFO("Filter mesh " << _mesh->getName() << " by bounding-box");
    Event e2("partition.filterMeshBB." + _mesh->getName());    
    mesh::Mesh filteredMesh("FilteredMesh", _dimensions, _mesh->isFlipNormals());
    filterMesh(filteredMesh, true);
    if ((_fromMapping.use_count() > 0 && _fromMapping->getOutputMesh()->vertices().size() > 0) ||
        (_toMapping.use_count() > 0 && _toMapping->getInputMesh()->vertices().size() > 0)) {
      // this rank has vertices at the coupling interface
      // then, also the filtered mesh should still have vertices
      std::string msg = "The re-partitioning completely filtered out the mesh " + _mesh->getName() + " received on this rank at the coupling interface. "
        "Most probably, the coupling interfaces of your coupled participants do not match geometry-wise. "
        "Please check your geometry setup again. Small overlaps or gaps are no problem. "
        "If your geometry setup is correct and if you have very different mesh resolutions on both sides, increasing the safety-factor "
        "of the decomposition strategy might be necessary.";
      CHECK(filteredMesh.vertices().size() > 0, msg);
    }

    DEBUG("Bounding box filter, filtered from " << _mesh->vertices().size() << " vertices to " << filteredMesh.vertices().size() << " vertices.");
    _mesh->clear();
    _mesh->addMesh(filteredMesh);
    _mesh->computeState();
    e2.stop();
  }
  

  // (2) Tag vertices 1st round (i.e. who could be owned by this rank)
  DEBUG("Tag vertices for filtering: 1st round.");
  // go to both meshes, vertex is tagged if already one mesh tags him
  if (_fromMapping.use_count() > 0)
    _fromMapping->tagMeshFirstRound();
  if (_toMapping.use_count() > 0)
    _toMapping->tagMeshFirstRound();

  // (3) Define which vertices are owned by this rank
  DEBUG("Create owner information.");
  createOwnerInformation();

  // (4) Tag vertices 2nd round (what should be filtered out)
  DEBUG("Tag vertices for filtering: 2nd round.");
  if (_fromMapping.use_count() > 0)
    _fromMapping->tagMeshSecondRound();
  if (_toMapping.use_count() > 0)
    _toMapping->tagMeshSecondRound();

  // (5) Filter mesh according to tag
  INFO("Filter mesh " << _mesh->getName() << " by mappings");
  Event e5("partition.filterMeshMappings" + _mesh->getName());
  mesh::Mesh filteredMesh("FilteredMesh", _dimensions, _mesh->isFlipNormals());
  filterMesh(filteredMesh, false);
  DEBUG("Mapping filter, filtered from " << _mesh->vertices().size() << " vertices to " << filteredMesh.vertices().size() << " vertices.");
  _mesh->clear();
  _mesh->addMesh(filteredMesh);
  _mesh->computeState();
  e5.stop();

  // (6) Compute distribution
  INFO("Feedback distribution for mesh " << _mesh->getName());
  Event e6("partition.feedbackMesh." + _mesh->getName());
  
  int              numberOfVertices = _mesh->vertices().size();
  std::vector<int> vertexIDs(numberOfVertices, -1);
  for (int i = 0; i < numberOfVertices; i++) {
    vertexIDs[i] = _mesh->vertices()[i].getGlobalIndex();
  }
  _mesh->getVertexDistribution()[utils::MasterSlave::_rank] = vertexIDs;

  remoteParComSize = vertexCounters.size();// number of other particpants ranks

  // This nested loop creats a fill in the localCommunicationbMap which shows this local rank needs which vertices from which rank of the other particpant
  for (auto &remoteVertex : vertexIDs) {
    for (int i=0; i <=remoteParComSize ; i++) {
      if (remoteVertex <= vertexCounters[i]) {
        localCommunicationMap[i].push_back(remoteVertex);
        i=remoteParComSize+1;
      }
    }
  }

  int globalNumberOfVertices = -1;
  _mesh->setGlobalNumberOfVertices(globalNumberOfVertices);
  _m2n->sendCommunicationMap(localCommunicationMap, *_mesh);  
  
  e6.stop();

  computeVertexOffsets();

}

bool ReceivedBoundingBox::CompareBoundingBox(mesh::Mesh::BoundingBox currentBB, mesh::Mesh::BoundingBox receivedBB)
{
  //int sizeofBB = currentBB.size();
  bool intersect=1;

  for (int i=0; i < _dimensions; i++) {

    if ((currentBB[i].first < receivedBB[i].first && currentBB[i].second < receivedBB[i].first) || (receivedBB[i].first < currentBB[i].first && receivedBB[i].second < currentBB[i].first) ) {

      intersect = 0;
      i=_dimensions;
    }
  }
  return intersect;
}

void ReceivedBoundingBox:: filterMesh(mesh::Mesh& filteredMesh, const bool filterByBB){
  TRACE(filterByBB);

  DEBUG("Bounding mesh. #vertices: " << _mesh->vertices().size()
               <<", #edges: " << _mesh->edges().size()
               <<", #triangles: " << _mesh->triangles().size() << ", rank: " << utils::MasterSlave::_rank);

  std::map<int, mesh::Vertex*> vertexMap;
  std::map<int, mesh::Edge*> edgeMap;
  int vertexCounter = 0;

  for (const mesh::Vertex& vertex : _mesh->vertices()) {

    if ((filterByBB && isVertexInBB(vertex)) || (not filterByBB && vertex.isTagged())){
      mesh::Vertex& v = filteredMesh.createVertex(vertex.getCoords());
      v.setGlobalIndex(vertex.getGlobalIndex());
      if(vertex.isTagged()) v.tag();
      v.setOwner(vertex.isOwner());
      vertexMap[vertex.getID()] = &v;
    }
    vertexCounter++;
  }

  // Add all edges formed by the contributing vertices
  for (mesh::Edge& edge : _mesh->edges()) {
    int vertexIndex1 = edge.vertex(0).getID();
    int vertexIndex2 = edge.vertex(1).getID();
    if (utils::contained(vertexIndex1, vertexMap) &&
        utils::contained(vertexIndex2, vertexMap)) {
      mesh::Edge& e = filteredMesh.createEdge(*vertexMap[vertexIndex1], *vertexMap[vertexIndex2]);
      edgeMap[edge.getID()] = &e;
    }
  }

  // Add all triangles formed by the contributing edges
  if (_dimensions==3) {
    for (mesh::Triangle& triangle : _mesh->triangles() ) {
      int edgeIndex1 = triangle.edge(0).getID();
      int edgeIndex2 = triangle.edge(1).getID();
      int edgeIndex3 = triangle.edge(2).getID();
      if (utils::contained(edgeIndex1, edgeMap) &&
          utils::contained(edgeIndex2, edgeMap) &&
          utils::contained(edgeIndex3, edgeMap)) {
        filteredMesh.createTriangle(*edgeMap[edgeIndex1],*edgeMap[edgeIndex2],*edgeMap[edgeIndex3]);
      }
    }
  }

  DEBUG("Filtered mesh. #vertices: " << filteredMesh.vertices().size()
               <<", #edges: " << filteredMesh.edges().size()
               <<", #triangles: " << filteredMesh.triangles().size() << ", rank: " << utils::MasterSlave::_rank);
}

void ReceivedBoundingBox::prepareBoundingBox(){
  TRACE(_safetyFactor);

  _bb.resize(_dimensions, std::make_pair(std::numeric_limits<double>::max(), std::numeric_limits<double>::lowest()));

  //create BB around both "other" meshes
  if (_fromMapping.use_count()>0) {
    auto other_bb = _fromMapping->getOutputMesh()->getBoundingBox();
    for (int d=0; d < _dimensions; d++) {
      if (_bb[d].first > other_bb[d].first) _bb[d].first = other_bb[d].first;
      if (_bb[d].second < other_bb[d].second) _bb[d].second = other_bb[d].second;
    }
  }
  if (_toMapping.use_count()>0) {
    auto other_bb = _toMapping->getInputMesh()->getBoundingBox();
    for (int d=0; d<_dimensions; d++) {
      if (_bb[d].first > other_bb[d].first) _bb[d].first = other_bb[d].first;
      if (_bb[d].second < other_bb[d].second) _bb[d].second = other_bb[d].second;
    }
  }

  //enlarge BB
  assertion(_safetyFactor>=0.0);

  double maxSideLength = 1e-6; // we need some minimum > 0 here

  for (int d=0; d<_dimensions; d++) {
    maxSideLength = std::max(maxSideLength, _bb[d].second - _bb[d].first);
  }
  for (int d=0; d<_dimensions; d++) {
    _bb[d].second += _safetyFactor * maxSideLength;
    _bb[d].first -= _safetyFactor * maxSideLength;
    DEBUG("Merged BoundingBox, dim: " << d << ", first: " << _bb[d].first << ", second: " << _bb[d].second);
  }


}

bool ReceivedBoundingBox::isVertexInBB(const mesh::Vertex& vertex){
  for (int d=0; d<_dimensions; d++) {
    if (vertex.getCoords()[d] < _bb[d].first || vertex.getCoords()[d] > _bb[d].second ) {
      return false;
    }
  }
  return true;
}

void ReceivedBoundingBox:: createOwnerInformation(){
  TRACE();

  if (utils::MasterSlave::_slaveMode) {
    int numberOfVertices = _mesh->vertices().size();
    utils::MasterSlave::_communication->send(numberOfVertices,0);

    if (numberOfVertices!=0) {
      std::vector<int> tags(numberOfVertices, -1);
      std::vector<int> globalIDs(numberOfVertices, -1);
      bool atInterface = false;
      for(int i=0; i<numberOfVertices; i++){
        globalIDs[i] = _mesh->vertices()[i].getGlobalIndex();
        if(_mesh->vertices()[i].isTagged()){
          tags[i] = 1;
          atInterface = true;
        }
        else{
          tags[i] = 0;
        }
      }
      DEBUG("My tags: " << tags);
      DEBUG("My global IDs: " << globalIDs);
      utils::MasterSlave::_communication->send(tags, 0);
      utils::MasterSlave::_communication->send(globalIDs, 0);
      utils::MasterSlave::_communication->send(atInterface, 0);

      std::vector<int> ownerVec(numberOfVertices, -1);
      utils::MasterSlave::_communication->receive(ownerVec, 0);
      DEBUG("My owner information: " << ownerVec);
      setOwnerInformation(ownerVec);
    }

  }


  else if (utils::MasterSlave::_masterMode) {
    std::vector<int> globalOwnerVec(_mesh->getGlobalNumberOfVertices(),0); //to temporary store which vertices already have an owner
    std::vector<std::vector<int> > slaveOwnerVecs; // the same per rank
    std::vector<std::vector<int> > slaveGlobalIDs; // global IDs per rank
    std::vector<std::vector<int> > slaveTags; // tag information per rank

    slaveOwnerVecs.resize(utils::MasterSlave::_size);
    slaveGlobalIDs.resize(utils::MasterSlave::_size);
    slaveTags.resize(utils::MasterSlave::_size);

    // fill master data

    bool masterAtInterface = false;
    slaveOwnerVecs[0].resize(_mesh->vertices().size());
    slaveGlobalIDs[0].resize(_mesh->vertices().size());
    slaveTags[0].resize(_mesh->vertices().size());
    for(size_t i=0; i<_mesh->vertices().size(); i++){
      slaveGlobalIDs[0][i] = _mesh->vertices()[i].getGlobalIndex();
      if(_mesh->vertices()[i].isTagged()){
        masterAtInterface = true;
        slaveTags[0][i] = 1;
      }
      else{
        slaveTags[0][i] = 0;
      }
    }
    DEBUG("My tags: " << slaveTags[0]);

    // receive slave data
    int ranksAtInterface = 0;
    if(masterAtInterface) ranksAtInterface++;

    for (int rank = 1; rank < utils::MasterSlave::_size; rank++){
      int localNumberOfVertices = -1;
      utils::MasterSlave::_communication->receive(localNumberOfVertices,rank);
    }
      
      DEBUG("Rank " << rank << " has " << localNumberOfVertices << " vertices.");
      slaveOwnerVecs[rank].resize(localNumberOfVertices, 0);
      slaveTags[rank].resize(localNumberOfVertices, -1);
      slaveGlobalIDs[rank].resize(localNumberOfVertices, -1);

      if (localNumberOfVertices!=0) {
        utils::MasterSlave::_communication->receive(slaveTags[rank], rank);
        utils::MasterSlave::_communication->receive(slaveGlobalIDs[rank], rank);
        DEBUG("Rank " << rank << " has this tags " << slaveTags[rank]);
        DEBUG("Rank " << rank << " has this global IDs " << slaveGlobalIDs[rank]);
        bool atInterface = false;
        utils::MasterSlave::_communication->receive(atInterface,rank);
        if(atInterface) ranksAtInterface++;
      }
    }

    // decide upon owners,
    int localGuess = _mesh->getGlobalNumberOfVertices() / ranksAtInterface; //guess for a decent load balancing
    //first round: every slave gets localGuess vertices
    for (int rank = 0; rank < utils::MasterSlave::_size; rank++){
      int counter = 0;
      for (size_t i=0; i < slaveOwnerVecs[rank].size(); i++) {
        if (globalOwnerVec[slaveGlobalIDs[rank][i]] == 0 && slaveTags[rank][i]==1) { // Vertex has no owner yet and rank could be owner
          slaveOwnerVecs[rank][i] = 1; // Now rank is owner
          globalOwnerVec[slaveGlobalIDs[rank][i]] = 1; //vertex now has owner
          counter++;
          if(counter==localGuess) break;
        }
      }
    }

    //second round: distribute all other vertices in a greedy way
    for (int rank = 0; rank < utils::MasterSlave::_size; rank++) {
      for(size_t i=0; i < slaveOwnerVecs[rank].size(); i++){
        if(globalOwnerVec[slaveGlobalIDs[rank][i]] == 0 && slaveTags[rank][i]==1){
          slaveOwnerVecs[rank][i] = 1;
          globalOwnerVec[slaveGlobalIDs[rank][i]] = rank + 1;
        }
      }
    }

    // send information back to slaves
    for (int rank = 1; rank < utils::MasterSlave::_size; rank++){
      int localNumberOfVertices = slaveTags[rank].size();
      if (localNumberOfVertices!=0) {
        utils::MasterSlave::_communication->send(slaveOwnerVecs[rank], rank);
      }
    }
    // master data
    DEBUG("My owner information: " << slaveOwnerVecs[0]);
    setOwnerInformation(slaveOwnerVecs[0]);
     

#     ifndef NDEBUG
    for(size_t i=0;i<globalOwnerVec.size();i++){
      if(globalOwnerVec[i]==0){
        WARN( "The Vertex with global index " << i << " of mesh: " << _mesh->getName()
                       << " was completely filtered out, since it has no influence on any mapping.");
          }
    }
#     endif

  }
}

void ReceivedBoundingBox:: setOwnerInformation(const std::vector<int> &ownerVec){
  size_t i = 0;
  for ( mesh::Vertex& vertex : _mesh->vertices() ){
    assertion(i<ownerVec.size());
    assertion(ownerVec[i]!=-1);
    vertex.setOwner(ownerVec[i]==1);
    i++;
  }
}


}}
