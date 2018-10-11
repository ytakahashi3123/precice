#include "partition/ProvidedBoundingBox.hpp"
#include "com/CommunicateBoundingBox.hpp"
#include "com/CommunicateMesh.hpp"
#include "com/Communication.hpp"
#include "utils/MasterSlave.hpp"
#include "m2n/M2N.hpp"
#include "utils/EventTimings.hpp"
#include "utils/Parallel.hpp"
#include "mesh/Mesh.hpp"
#include "mapping/Mapping.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/Vertex.hpp"
#include "mesh/Edge.hpp"
#include "mesh/Triangle.hpp"


using precice::utils::Event;

namespace precice {
namespace partition {

logging::Logger ProvidedBoundingBox:: _log ( "precice::partition::ProvidedBoundingBox" );

ProvidedBoundingBox::ProvidedBoundingBox
(mesh::PtrMesh mesh,
 bool hasToSend,
 double safetyFactor)
:
    Partition (mesh),
    _hasToSend(hasToSend),
    _bb(mesh->getBoundingBox()),
    _dimensions(mesh->getDimensions()),
    _safetyFactor(safetyFactor)
{}

void ProvidedBoundingBox::communicateBoundingBox()
{

  if (_hasToSend) { 
  
    Event e1("creat and gather bounding box");

    // each rank sends its bb to master
    if (utils::MasterSlave::_slaveMode) {//slave
      com::CommunicateBoundingBox(utils::MasterSlave::_communication).sendBoundingBox(_bb, 0); 
    }
    else{ // Master

      assertion(utils::MasterSlave::_rank==0);
      assertion(utils::MasterSlave::_size>1);

      // master stores its bb into gloabalbb
      _globalBB[0] = _bb;

      // master receives bbs from slaves and store them in global bb
      if (utils::MasterSlave::_size>1) {  
      for (int rankSlave = 1; rankSlave < utils::MasterSlave::_size; rankSlave++) {
        com::CommunicateBoundingBox(utils::MasterSlave::_communication).receiveBoundingBox(_bb, rankSlave);
        
        DEBUG("From slave " << rankSlave << ", bounding mesh: " << _bb[0].first
                     << ", " << _bb[0].second << " and " << _bb[1].first << ", " << _bb[1].second);
        
        _globalBB[rankSlave] = _bb;
      }
      }
    }
    
    e1.stop();

    //master sends global bb to the other master
    Event e2("send global Bounding Box");
    if (utils::MasterSlave::_masterMode) {
      _m2n->getMasterCommunication()->send(utils::MasterSlave::_size , 0);
      com::CommunicateBoundingBox(_m2n->getMasterCommunication()).sendBoundingBoxMap(_globalBB,0);
    }
    e2.stop();
  }

   
}
void ProvidedBoundingBox::computeBoundingBox()
{
  int numberOfVertices = _mesh->vertices().size();

    if (not utils::MasterSlave::_slaveMode) {//Master

      assertion(utils::MasterSlave::_size>1);
      int vertexCounter = 0;

      //master receives other partition communicator size and also a feedback which is a map:  list of other participant ranks -> connected ranks at this participant  
      _m2n->getMasterCommunication()->receive(remoteParComSize, 0);
      utils::MasterSlave::_communication->broadcast(remoteParComSize);
      for (int i=0; i < remoteParComSize; i++) {
        std::vector<int> test;
        test.push_back(-1);
        received_feedbackMap[i]=test;
      }
      com::CommunicateBoundingBox(_m2n->getMasterCommunication()).receiveFeedbackMap(received_feedbackMap, 0 ); 
      com::CommunicateBoundingBox(utils::MasterSlave::_communication).broadcastSendFeedbackMap(received_feedbackMap);     


      // master checks which ranks are connected to it
    for (auto &other_rank : received_feedbackMap) {
        for (auto &included_ranks: other_rank.second) {
          if (utils::Parallel::getProcessRank() == included_ranks) {
            connected_ranks.push_back(other_rank.first);
            _mesh->getCommunicationMap()[other_rank.first].push_back(-1);
          }
        }
      }
    
    //add master vertices to vertex counter which is a vector that shows each ranks maximum vertex id
    for(int i=0; i<numberOfVertices; i++){
      _mesh->getVertexDistribution()[0].push_back(vertexCounter);
      _mesh->vertices()[i].setGlobalIndex(vertexCounter);
      vertexCounter++;
    }
    
    vertexCounters.push_back(vertexCounter);
    

    for (int rankSlave = 1; rankSlave < utils::MasterSlave::_size; rankSlave++){
      utils::MasterSlave::_communication->receive(numberOfVertices,rankSlave);
      utils::MasterSlave::_communication->send(vertexCounter,rankSlave);      

      for(int i=0; i<numberOfVertices; i++){
        _mesh->getVertexDistribution()[rankSlave].push_back(vertexCounter);
        vertexCounter++;
      }
      vertexCounters[rankSlave]=vertexCounter;

    }
   
    _mesh->setGlobalNumberOfVertices(vertexCounter);

    for (int j = 0; j < utils::MasterSlave::_size; j++) {      
    utils::MasterSlave::_communication->broadcast(vertexCounters[j]);
    }
    
    }
    else{ // Slave

      utils::MasterSlave::_communication->broadcast(remoteParComSize, 0);
      for (int i=0; i < remoteParComSize; i++) {
        std::vector<int> test;
        test.push_back(i);
        received_feedbackMap[i]=test;
      }
      com::CommunicateBoundingBox(utils::MasterSlave::_communication).broadcastReceiveFeedbackMap(received_feedbackMap);

      for (auto &other_rank : received_feedbackMap) {
        for (auto &included_ranks: other_rank.second) {
          if (utils::Parallel::getProcessRank() == included_ranks) {
            connected_ranks.push_back(other_rank.first);
            _mesh->getCommunicationMap()[other_rank.first].push_back(-1);
          }
        }
      }
      
    int globalVertexCounter = 1;
    utils::MasterSlave::_communication->send(numberOfVertices,0);
    utils::MasterSlave::_communication->receive(globalVertexCounter,0);
    for(int i=0; i<numberOfVertices; i++){
      _mesh->vertices()[i].setGlobalIndex(globalVertexCounter+i);
    }
    int globalNumberOfVertices = 1;

    vertexCounters.resize(utils::MasterSlave::_size);
    for (int j = 0; j < utils::MasterSlave::_size; j++) {      
      utils::MasterSlave::_communication->broadcast(vertexCounters[j], 0);
    }
    
    globalNumberOfVertices = vertexCounters.back();
    assertion(globalNumberOfVertices!=-1);
    _mesh->setGlobalNumberOfVertices(globalNumberOfVertices);
    
    }
          
  createOwnerInformation();

  computeVertexOffsets();

  if (utils::MasterSlave::_masterMode)
  {
    for (int i = 0; i <utils::MasterSlave::_size ; i++) {
      _m2n->getMasterCommunication()->send(vertexCounters[i], 0);
    }
  }
}

void ProvidedBoundingBox::communicate()
{

  // each rank send its mesh partition to connected ranks at other participant
  _m2n->broadcastSendLocalMesh(*_mesh);
}

void ProvidedBoundingBox::compute()
{
  // each rank receives its final communication map
//  _m2n->receiveCommunicationMap(_mesh->getCommunicationMap(), *_mesh);
}

void ProvidedBoundingBox::createOwnerInformation()
{
  TRACE();
  for(mesh::Vertex& v : _mesh->vertices()){
    v.setOwner(true);
  }
}

}}
