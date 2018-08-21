#include "partition/ReceivedBoundingBox.hpp"
#include "com/CommunicateBoundingBox.hpp"
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
#include "utils/Globals.hpp"
#include "utils/MasterSlave.hpp"
#include <vector>
#include <map>

using precice::utils::Event;

namespace precice {
namespace partition {

logging::Logger ReceivedBoundingBox:: _log ( "precice::partition::ReceivedBoundingBox" );

ReceivedBoundingBox::ReceivedBoundingBox
(
  mesh::PtrMesh mesh, double safetyFactor, mesh::Mesh::BoundingBoxMap globalBB)
:
  Partition (mesh),
  /* _bb(mesh->getDimensions(), std::make_pair(std::numeric_limits<double>::max(), std::numeric_limits<double>::lowest())),*/
  _bb(mesh->getBoundingBox()),
  _dimensions(mesh->getDimensions()),
  _safetyFactor(safetyFactor),
  _globalBB(globalBB)
{}

void ReceivedBoundingBox::communicate()
{
  TRACE();
  Event e("receive global bounding box");  
  if (not utils::MasterSlave::_slaveMode) {     
  com::CommunicateBoundingBox(_m2n->getMasterCommunication()).receiveBoundingBoxMap(_globalBB, 0 );
  }
}

void ReceivedBoundingBox::compute()
{
  // handle coupling mode first (i.e. serial participant)
  if (not utils::MasterSlave::_slaveMode && not utils::MasterSlave::_masterMode)
  { //coupling mode
  //some code here
  }

  Event e1("broadcast received bounding box");
   

    if (utils::MasterSlave::_slaveMode)
    {
      mesh::Mesh::BoundingBoxMap received_globalBB;
      com::CommunicateBoundingBox(utils::MasterSlave::_communication).broadcastReceiveBoundingBoxMap(received_globalBB);
      numberOfVertices = _mesh->vertices().size();
      prepareBoundingBox();
      feedback.push_back(-1);

      if (numberOfVertices>0) {
        for (auto &other_rank: received_globalBB)
        {
          if (CompareBoundingBox(_bb,other_rank.second)) {
            feedback.push_back(other_rank.first);
          }
        }
      }
 

     //send feedback to master
     utils::MasterSlave::_communication->send(feedback, 0);
     
    }
    else 
    { // Master

      assertion(utils::MasterSlave::_rank==0);
      assertion(utils::MasterSlave::_size>1);       
      com::CommunicateBoundingBox(utils::MasterSlave::_communication).broadcastSendBoundingBoxMap(_globalBB);
      numberOfVertices = _mesh->vertices().size();

      prepareBoundingBox();
      feedbackMap[0].push_back(-1);
      if (numberOfVertices>0) {           
        for (auto &other_rank: _globalBB)
        {
          if (CompareBoundingBox(_bb,other_rank.second)) {
            feedbackMap[0].push_back(other_rank.first);
          }
        }
      }
      
      for (int rank_slave=1; rank_slave < utils::MasterSlave::_size ; rank_slave++) {     
        utils::MasterSlave::_communication->receive(feedback, rank_slave);
        feedbackMap[rank_slave]=feedback;        
      }

      com::CommunicateBoundingBox(_m2n->getMasterCommunication()).sendFeedbackMap(feedbackMap, 0 ); // @Amin: create this method!

      
    }
    e1.stop();
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

void ReceivedBoundingBox::prepareBoundingBox(){
  TRACE(_safetyFactor);

  _bb.resize(_dimensions, std::make_pair(std::numeric_limits<double>::max(), std::numeric_limits<double>::lowest()));

  //enlarge BB
  assertion(_safetyFactor>=0.0);

  double maxSideLength = 1e-6; // we need some minimum > 0 here

  for (int d=0; d<_dimensions; d++) {
    maxSideLength = std::max(maxSideLength, _bb[d].second - _bb[d].first);
  }
  for (int d=0; d<_dimensions; d++) {
    _bb[d].second += _safetyFactor * maxSideLength;
    _bb[d].first -= _safetyFactor * maxSideLength;    
  }
}


void ReceivedBoundingBox::createOwnerInformation()
{}

}}
