#ifndef PRECICE_NO_MPI
#include "testing/Testing.hpp"
#include "testing/Fixtures.hpp"

#include "partition/ProvidedBoundingBox.hpp"
#include "partition/ReceivedBoundingBox.hpp"

#include "utils/Parallel.hpp"
#include "com/MPIDirectCommunication.hpp"
#include "com/SocketCommunication.hpp"
#include "com/CommunicateBoundingBox.hpp"
#include "m2n/M2N.hpp"
#include "utils/MasterSlave.hpp"
#include "m2n/GatherScatterComFactory.hpp"
#include "m2n/PointToPointCommunication.hpp"
#include "mapping/NearestNeighborMapping.hpp"
#include "mapping/NearestProjectionMapping.hpp"
#include "mapping/PetRadialBasisFctMapping.hpp"
#include "mesh/Vertex.hpp"
#include "mesh/Edge.hpp"
#include "mesh/Triangle.hpp"
#include "mesh/Quad.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/PropertyContainer.hpp"
#include "mesh/Data.hpp"

using namespace precice;
using namespace partition;

BOOST_AUTO_TEST_SUITE(PartitionTests)
BOOST_AUTO_TEST_SUITE(ReceivedBoundingBoxTests)

void setupParallelEnvironment(m2n::PtrM2N m2n) {
  assertion(utils::Parallel::getCommunicatorSize() == 4);

  com::PtrCommunication masterSlaveCom =
        com::PtrCommunication(new com::MPIDirectCommunication());
  utils::MasterSlave::_communication = masterSlaveCom;

  utils::Parallel::synchronizeProcesses();

  if (utils::Parallel::getProcessRank() == 0){ //NASTIN
    utils::Parallel::splitCommunicator( "Fluid" );
    m2n->acceptMasterConnection ( "Fluid", "SolidMaster");
    utils::MasterSlave::_slaveMode = false;
    utils::MasterSlave::_masterMode = false;
  }
  else if(utils::Parallel::getProcessRank() == 1){//Master
    utils::Parallel::splitCommunicator( "SolidMaster" );
    m2n->requestMasterConnection ( "Fluid", "SolidMaster" );
    utils::MasterSlave::_rank = 0;
    utils::MasterSlave::_size = 3;
    utils::MasterSlave::_slaveMode = false;
    utils::MasterSlave::_masterMode = true;
  }
  else if(utils::Parallel::getProcessRank() == 2){//Slave1
    utils::Parallel::splitCommunicator( "SolidSlaves");
    utils::MasterSlave::_rank = 1;
    utils::MasterSlave::_size = 3;
    utils::MasterSlave::_slaveMode = true;
    utils::MasterSlave::_masterMode = false;
  }
  else if(utils::Parallel::getProcessRank() == 3){//Slave2
    utils::Parallel::splitCommunicator( "SolidSlaves");
    utils::MasterSlave::_rank = 2;
    utils::MasterSlave::_size = 3;
    utils::MasterSlave::_slaveMode = true;
    utils::MasterSlave::_masterMode = false;
  }

  if(utils::Parallel::getProcessRank() == 1){//Master
    masterSlaveCom->acceptConnection ( "SolidMaster", "SolidSlaves");
    masterSlaveCom->setRankOffset(1);
  }
  else if(utils::Parallel::getProcessRank() == 2){//Slave1
    masterSlaveCom->requestConnection( "SolidMaster", "SolidSlaves", 0, 2 );
  }
  else if(utils::Parallel::getProcessRank() == 3){//Slave2
    masterSlaveCom->requestConnection( "SolidMaster", "SolidSlaves", 1, 2 );
  }
  }



void tearDownParallelEnvironment(){
  utils::MasterSlave::_communication = nullptr;
  utils::MasterSlave::reset();
  utils::Parallel::synchronizeProcesses();
  utils::Parallel::clearGroups();
  mesh::Mesh::resetGeometryIDsGlobally();
  mesh::Data::resetDataCount();
  utils::Parallel::setGlobalCommunicator(utils::Parallel::getCommunicatorWorld());
}

void createNastinMesh2D(mesh::PtrMesh pNastinMesh){
  int dimensions = 2;
  assertion(pNastinMesh.use_count()>0);
  assertion(pNastinMesh->getDimensions()==dimensions);
  Eigen::VectorXd position(dimensions);

  position << -3.0,-10.0;
  mesh::Vertex& v1 = pNastinMesh->createVertex(position);
  v1.setGlobalIndex(0);
  position << 10.0, 1.95;
  mesh::Vertex& v2 = pNastinMesh->createVertex(position);
  v2.setGlobalIndex(1);
  position << 0.0, 2.1;
  mesh::Vertex& v3 = pNastinMesh->createVertex(position);
  v3.setGlobalIndex(2);
  position << 0.0, 4.5;
  mesh::Vertex& v4 = pNastinMesh->createVertex(position);
  v4.setGlobalIndex(3);
  position << 0.0, 5.95;
  mesh::Vertex& v5 = pNastinMesh->createVertex(position);
  v5.setGlobalIndex(4);
  position << 0.0, 6.1;
  mesh::Vertex& v6 = pNastinMesh->createVertex(position);
  v6.setGlobalIndex(5);
  pNastinMesh->createEdge(v1,v2);
  pNastinMesh->createEdge(v2,v3);
  pNastinMesh->createEdge(v3,v4);
  pNastinMesh->createEdge(v4,v5);
  pNastinMesh->createEdge(v5,v6);
}

void createSolidzMesh2D(mesh::PtrMesh pSolidzMesh){
  int dimensions = 2;
  assertion(pSolidzMesh.use_count()>0);
  assertion(pSolidzMesh->getDimensions()==dimensions);

  if(utils::Parallel::getProcessRank() == 1){

    Eigen::VectorXd position(dimensions);
    position << 0.0, 0.0;
    pSolidzMesh->createVertex(position);
    position << 0.0, 2.0;
    pSolidzMesh->createVertex(position);
    position << 0.1, 3.0;
    pSolidzMesh->createVertex(position);
  }
  else if(utils::Parallel::getProcessRank() == 2){
    // not at interface
  }
  else if(utils::Parallel::getProcessRank() == 3){

    Eigen::VectorXd position(dimensions);
    position << -1.0, 4.0;
    pSolidzMesh->createVertex(position);
    position << 2.0, 3.0;
    pSolidzMesh->createVertex(position);
  }
}

void createNastinMesh3D(mesh::PtrMesh pNastinMesh){
  int dimensions = 3;
  Eigen::VectorXd position(dimensions);
  assertion(pNastinMesh.use_count()>0);
  assertion(pNastinMesh->getDimensions()==dimensions);
  
  position << 0.0, 0.0, -0.1;
  mesh::Vertex& v1 = pNastinMesh->createVertex(position);
  v1.setGlobalIndex(0);
  position << -1.0, 0.0, 0.0;
  mesh::Vertex& v2 = pNastinMesh->createVertex(position);
  v2.setGlobalIndex(1);
  position << 1.0, 0.0, 0.0;
  mesh::Vertex& v3 = pNastinMesh->createVertex(position);
  v3.setGlobalIndex(2);
  position << 0.0, -1.0, 0.0;
  mesh::Vertex& v4 = pNastinMesh->createVertex(position);
  v4.setGlobalIndex(3);
  position << 0.0, 1.0, 0.0;
  mesh::Vertex& v5 =pNastinMesh->createVertex(position);
  v5.setGlobalIndex(4);
  mesh::Edge& e1 = pNastinMesh->createEdge(v1,v2);
  mesh::Edge& e2 = pNastinMesh->createEdge(v2,v4);
  mesh::Edge& e3 = pNastinMesh->createEdge(v4,v1);
  mesh::Edge& e4 = pNastinMesh->createEdge(v1,v3);
  mesh::Edge& e5 = pNastinMesh->createEdge(v3,v5);
  mesh::Edge& e6 = pNastinMesh->createEdge(v5,v1);
  pNastinMesh->createTriangle(e1,e2,e3);
  pNastinMesh->createTriangle(e4,e5,e6);
}



void createSolidzMesh3D(mesh::PtrMesh pSolidzMesh){
  int dimensions = 3;
  assertion(pSolidzMesh.use_count()>0);
  assertion(pSolidzMesh->getDimensions()==dimensions);

  if(utils::Parallel::getProcessRank() == 1){//Master

    Eigen::VectorXd position(dimensions);
    position << -1.0, -1.0, 0.0;
    pSolidzMesh->createVertex(position);
    position << -0.75, -0.75, 0.5;
    pSolidzMesh->createVertex(position);
  }
  else if(utils::Parallel::getProcessRank() == 2){//Slave1
    // slave1 not at interface
  }
  else if(utils::Parallel::getProcessRank() == 3){//Slave2

    Eigen::VectorXd position(dimensions);
    position << 0.0, 0.0, -1.0;
    pSolidzMesh->createVertex(position);
    position << 0.5, 0.5, 0.0;
    pSolidzMesh->createVertex(position);
  }
}


BOOST_AUTO_TEST_CASE(TestReceivedBoundingBox2D, * testing::OnSize(4))
{
  com::PtrCommunication participantCom =
      com::PtrCommunication(new com::SocketCommunication());
  m2n::DistributedComFactory::SharedPointer distrFactory = m2n::DistributedComFactory::SharedPointer(
      new m2n::GatherScatterComFactory(participantCom));
  m2n::PtrM2N m2n = m2n::PtrM2N(new m2n::M2N(participantCom, distrFactory));

  setupParallelEnvironment(m2n);

  int dimensions = 2;
  bool flipNormals = false;
  Eigen::VectorXd offset = Eigen::VectorXd::Zero(dimensions);

  if (utils::Parallel::getProcessRank() == 0){ //NASTIN
    mesh::PtrMesh pNastinMesh(new mesh::Mesh("NastinMesh", dimensions, flipNormals));    
    double safetyFactor = 0;
    int remote_size;
    Eigen::VectorXd position(dimensions);
    createNastinMesh2D(pNastinMesh);
    pNastinMesh->computeState();   
    mesh::Mesh::BoundingBoxMap send_globalBB;    
    mesh::Mesh::BoundingBox localBB;
    mesh::Mesh::FeedbackMap feedbackMap;

    //prepare demo golabl boundingbox to send to receiving partition
    for (int i=0; i < 3; i++) {
      for (int j=0; j < 2; j++) {   
      localBB.push_back(std::make_pair(i,i+1));
      }
      send_globalBB[i]=localBB;
      localBB.clear();
    }


    for (int i=0; i < 3; i++) {
        std::vector<int> test;
        test.push_back(-1);
        feedbackMap[i]=test;
    }

    std::vector<int> vertexcounters;
    vertexcounters.push_back(6);
    int com_size = 3;
    m2n->getMasterCommunication()->send(com_size , 0);
    com::CommunicateBoundingBox(m2n->getMasterCommunication()).sendBoundingBoxMap(send_globalBB,0);
    m2n->getMasterCommunication()->receive(remote_size, 0);
    com::CommunicateBoundingBox(m2n->getMasterCommunication()).receiveFeedbackMap(feedbackMap,0);
    for (int i=0; i < com_size; i++) {
      m2n->getMasterCommunication()->send(i+2, 0);
    }
  }

  else{//SOLIDZ
    mesh::PtrMesh pSolidzMesh(new mesh::Mesh("SolidzMesh", dimensions, flipNormals));
    mesh::PtrMesh pSolidzMesh_received(new mesh::Mesh("SolidzMesh", dimensions, flipNormals));

    mapping::PtrMapping boundingFromMapping = mapping::PtrMapping(new mapping::NearestNeighborMapping(mapping::Mapping::CONSISTENT, dimensions));
    mapping::PtrMapping boundingToMapping = mapping::PtrMapping(new mapping::NearestNeighborMapping(mapping::Mapping::CONSERVATIVE, dimensions));
    boundingFromMapping->setMeshes(pSolidzMesh_received, pSolidzMesh);
    boundingToMapping->setMeshes(pSolidzMesh, pSolidzMesh_received);


    if(utils::Parallel::getProcessRank() == 1){//Master
      Eigen::VectorXd position(dimensions);
      position << 2.25, 2.61;
      pSolidzMesh->createVertex(position);
      position << 2.67, 2.95;
      pSolidzMesh->createVertex(position);
      position << 2.54, 2.10;
      pSolidzMesh->createVertex(position);
    }   
  
    else if(utils::Parallel::getProcessRank() == 2){//Slave1
      Eigen::VectorXd position(dimensions);
      position << 1.3, 1.5;
      mesh::Vertex& v3 = pSolidzMesh->createVertex(position);
      position << 1.65, 1.52;
      mesh::Vertex& v4 = pSolidzMesh->createVertex(position);
      pSolidzMesh->createEdge(v3,v4);
    }
    else if(utils::Parallel::getProcessRank() == 3){//Slave2
      Eigen::VectorXd position(dimensions);
      position << 0.2, 0.7;
      mesh::Vertex& v5 = pSolidzMesh->createVertex(position);
      position << 0.8, 0.4;
      mesh::Vertex& v6 = pSolidzMesh->createVertex(position); 
      pSolidzMesh->createEdge(v5,v6);
    }
 
    pSolidzMesh->computeState();
    
   
    double safetyFactor = 0;    
    mesh::Mesh::BoundingBoxMap globalBBtoCheck;
    mesh::Mesh::BoundingBox localBB;
      
    ReceivedBoundingBox part(pSolidzMesh_received, safetyFactor, ReceivedBoundingBox::FILTER_FIRST);
    part.setM2N(m2n);
    part.setFromMapping(boundingFromMapping);
    part.setToMapping(boundingToMapping);
    part.communicateBoundingBox();
    part.computeBoundingBox();
    
    for (int i=0; i < 3; i++) {    
      for (int j = 0; j < 2; j++) { 
        localBB.push_back(std::make_pair(i,i+1));
      }
      globalBBtoCheck[i]=localBB;
      localBB.clear();
    }

    // if (utils::Parallel::getProcessRank() == 1){
       BOOST_TEST(part.remoteParComSize==3);
    // }
    for (int i=0; i<3 ; i++) {                
      BOOST_TEST(part._globalBB[i] == globalBBtoCheck[i]);
    }
    
    if(utils::Parallel::getProcessRank()==1){      
      BOOST_TEST(part.feedbackMap[0][0]==2);
      BOOST_TEST(part.feedbackMap[1][0]==1);
      BOOST_TEST(part.feedbackMap[2][0]==0);
    } else if(utils::Parallel::getProcessRank()==2){      
      BOOST_TEST(part.feedback.size()==1);
      BOOST_TEST(part.feedback[0]==1);
    } else if(utils::Parallel::getProcessRank()==3){      
      BOOST_TEST(part.feedback.size()==1);
      BOOST_TEST(part.feedback[0]==0);
    }
   

  }
  
  tearDownParallelEnvironment();
}


BOOST_AUTO_TEST_CASE(TestReceivedBoundingBox3D, * testing::OnSize(4))
{
  com::PtrCommunication participantCom =
      com::PtrCommunication(new com::MPIDirectCommunication());
  m2n::DistributedComFactory::SharedPointer distrFactory = m2n::DistributedComFactory::SharedPointer(
      new m2n::GatherScatterComFactory(participantCom));
  m2n::PtrM2N m2n = m2n::PtrM2N(new m2n::M2N(participantCom, distrFactory));

  setupParallelEnvironment(m2n);

  int dimensions = 3;
  bool flipNormals = false;
  int remote_size=0;
  Eigen::VectorXd offset = Eigen::VectorXd::Zero(dimensions);

  if (utils::Parallel::getProcessRank() == 0){ //NASTIN

    mesh::Mesh::BoundingBoxMap send_globalBB;    
    mesh::Mesh::BoundingBox localBB;
    mesh::Mesh::FeedbackMap feedbackMap;

    //prepare demo golabl boundingbox to send to receiving partition
    for (int i=0; i < 3; i++) {
      for (int j=0; j < dimensions; j++) {
        localBB.push_back(std::make_pair(i,i+1));
      }      
      send_globalBB[i]=localBB;
      localBB.clear();
    }

    int com_size = 3;
    m2n->getMasterCommunication()-> send(com_size , 0);
    com::CommunicateBoundingBox(m2n->getMasterCommunication()).sendBoundingBoxMap(send_globalBB,0);     
    m2n->getMasterCommunication()->receive(remote_size, 0);
    com::CommunicateBoundingBox(m2n->getMasterCommunication()).receiveFeedbackMap(feedbackMap,0);
    for (int i=0; i < com_size; i++) {
      m2n->getMasterCommunication()->send(i+2, 0);
    }
  }

  else{//SOLIDZ
    mesh::PtrMesh pSolidzMesh(new mesh::Mesh("SolidzMesh", dimensions, flipNormals));
    mesh::PtrMesh pSolidzMesh_received(new mesh::Mesh("SolidzMesh", dimensions, flipNormals));

    mapping::PtrMapping boundingFromMapping = mapping::PtrMapping(new mapping::NearestNeighborMapping(mapping::Mapping::CONSISTENT, dimensions));
    mapping::PtrMapping boundingToMapping = mapping::PtrMapping(new mapping::NearestNeighborMapping(mapping::Mapping::CONSERVATIVE, dimensions));
    boundingFromMapping->setMeshes(pSolidzMesh_received, pSolidzMesh);
    boundingToMapping->setMeshes(pSolidzMesh, pSolidzMesh_received);

    if(utils::Parallel::getProcessRank() == 1){//Master
      Eigen::VectorXd position(dimensions);
      position << 2.25, 2.61, 2.24;
      pSolidzMesh->createVertex(position);
      position << 2.67, 2.95, 2.65;
      pSolidzMesh->createVertex(position);
      position << 2.54, 2.10, 2.32;
      pSolidzMesh->createVertex(position);
    }   
  
    else if(utils::Parallel::getProcessRank() == 2){//Slave1
      Eigen::VectorXd position(dimensions);
      position << 1.3, 1.5, 1.7;
      mesh::Vertex& v3 = pSolidzMesh->createVertex(position);
      position << 1.65, 1.52, 1.95;
      mesh::Vertex& v4 = pSolidzMesh->createVertex(position);
      pSolidzMesh->createEdge(v3,v4);
    }
    else if(utils::Parallel::getProcessRank() == 3){//Slave2
      Eigen::VectorXd position(dimensions);
      position << 0.2, 0.7, 0.8;
      mesh::Vertex& v5 = pSolidzMesh->createVertex(position);
      position << 0.8, 0.4, 0.05;
      mesh::Vertex& v6 = pSolidzMesh->createVertex(position); 
      pSolidzMesh->createEdge(v5,v6);
    }
 
    pSolidzMesh->computeState();
    
   
    double safetyFactor = 0;    
    mesh::Mesh::BoundingBoxMap globalBBtoCheck;
    mesh::Mesh::BoundingBox localBB;
      
    ReceivedBoundingBox part(pSolidzMesh_received, safetyFactor, ReceivedBoundingBox::FILTER_FIRST);
    part.setM2N(m2n);
    part.setFromMapping(boundingFromMapping);
    part.setToMapping(boundingToMapping);
    part.communicateBoundingBox();
    part.computeBoundingBox();
    
    for (int i=0; i < 3; i++) {
      for (int j=0; j < dimensions; j++) {
        localBB.push_back(std::make_pair(i,i+1));
      }      
      globalBBtoCheck[i]=localBB;
      localBB.clear();
    }

    BOOST_TEST(part.remoteParComSize==3);
    for (int i=0; i<3 ; i++) {                
      BOOST_TEST(part._globalBB[i] == globalBBtoCheck[i]);
    }

    if(utils::Parallel::getProcessRank()==1){      
      BOOST_TEST(part._bb[0].first==2.25);
      BOOST_TEST(part._bb[0].second==2.67);
      BOOST_TEST(part._bb[1].first==2.1);
      BOOST_TEST(part._bb[1].second==2.95);
      BOOST_TEST(part._bb[2].first==2.24);
      BOOST_TEST(part._bb[2].second==2.65);
    }else if(utils::Parallel::getProcessRank()==2){      
      BOOST_TEST(part._bb[0].first==1.3);
      BOOST_TEST(part._bb[0].second==1.65);
      BOOST_TEST(part._bb[1].first==1.5);
      BOOST_TEST(part._bb[1].second==1.52);
      BOOST_TEST(part._bb[2].first==1.7);
      BOOST_TEST(part._bb[2].second==1.95);
    }else if(utils::Parallel::getProcessRank()==3){      
      BOOST_TEST(part._bb[0].first==0.2);
      BOOST_TEST(part._bb[0].second==0.8);
      BOOST_TEST(part._bb[1].first==0.4);
      BOOST_TEST(part._bb[1].second==0.7);
      BOOST_TEST(part._bb[2].first==0.05);
      BOOST_TEST(part._bb[2].second==0.8);
    }
    

    if(utils::Parallel::getProcessRank()==1){      
      BOOST_TEST(part.feedbackMap.size()==3);
      BOOST_TEST(part.feedbackMap[0][0]==2);
      BOOST_TEST(part.feedbackMap[1][0]==1);
      BOOST_TEST(part.feedbackMap[2][0]==0);
      BOOST_TEST(part.feedback.size()==1);
    } else if(utils::Parallel::getProcessRank()==2){
      BOOST_TEST(part.feedback.size()==1);
      BOOST_TEST(part.feedback[0]==1);
    } else if(utils::Parallel::getProcessRank()==3){      
      BOOST_TEST(part.feedback.size()==1);
      BOOST_TEST(part.feedback[0]==0);
    }

    }
  tearDownParallelEnvironment();
}


BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()
#endif // PRECICE_NO_MPI
