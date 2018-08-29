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

void ProvidedBoundingBox::communicate()
{

  if (_hasToSend) { 
  
    Event e1("creat and gather bounding box");

    if (utils::MasterSlave::_slaveMode) {//slave
      com::CommunicateBoundingBox(utils::MasterSlave::_communication).sendBoundingBox(_bb, 0); 
    }
    else{ // Master

      assertion(utils::MasterSlave::_rank==0);
      assertion(utils::MasterSlave::_size>1);

      _globalBB[0] = _bb;
       
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

    Event e2("send global Bounding Box");
    if (utils::MasterSlave::_masterMode) {
      _m2n->getMasterCommunication()->send(utils::MasterSlave::_size , 0);
      com::CommunicateBoundingBox(_m2n->getMasterCommunication()).sendBoundingBoxMap(_globalBB,0);
    }
    e2.stop();
  }

   
}
void ProvidedBoundingBox::compute()
{
  int numberOfVertices = _mesh->vertices().size();

    if (not utils::MasterSlave::_slaveMode) {//Master

      assertion(utils::MasterSlave::_size>1);
      int vertexCounter = 0;

      _m2n->getMasterCommunication()->receive(remoteParComSize, 0);
      utils::MasterSlave::_communication->broadcast(remoteParComSize);
      for (int i=0; i < remoteParComSize; i++) {
        std::vector<int> test;
        test.push_back(i);
        received_feedbackMap[i]=test;
      }
 
      com::CommunicateBoundingBox(_m2n->getMasterCommunication()).receiveFeedbackMap(received_feedbackMap, 0 ); 
      com::CommunicateBoundingBox(utils::MasterSlave::_communication).broadcastSendFeedbackMap(received_feedbackMap);     

    for (auto &other_rank : received_feedbackMap) {
        for (auto &included_ranks: other_rank.second) {
          if (utils::Parallel::getProcessRank() == included_ranks) {
            connected_ranks.push_back(other_rank.first);
            _mesh->getCommunicationMap()[other_rank.first].push_back(-1);
          }
        }
      }
    
    //add master vertices
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
}

void ProvidedBoundingBox::communicatePartition()
{
  if (utils::MasterSlave::_masterMode)
  {
   _m2n->getMasterCommunication()->send(vertexCounters, 0);  
  }
  _m2n->sendMesh(*_mesh);
}

void ProvidedBoundingBox::computePartition()
{
  _m2n->receiveCommunicationMap(_mesh->getCommunicationMap(), *_mesh);
}

void ProvidedBoundingBox::createOwnerInformation()
{
  TRACE();
  for(mesh::Vertex& v : _mesh->vertices()){
    v.setOwner(true);
  }
}

}}
