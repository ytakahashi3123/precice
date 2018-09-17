#ifndef PRECICE_NO_MPI
#include "testing/Testing.hpp"
#include "testing/Fixtures.hpp"

#include "partition/ProvidedBoundingBox.hpp"
#include "partition/ReceivedBoundingBox.hpp"
#include "partition/SharedPointer.hpp"
#include "utils/Parallel.hpp"
#include "com/MPIDirectCommunication.hpp"
#include "com/MPIPortsCommunication.hpp"
#include "com/MPIPortsCommunicationFactory.hpp"
#include "com/SocketCommunication.hpp"
#include "com/SocketCommunicationFactory.hpp"
#include "com/CommunicateBoundingBox.hpp"
#include "m2n/PointToPointCommunication.hpp"
#include "m2n/M2N.hpp"
#include "utils/MasterSlave.hpp"
#include "m2n/GatherScatterComFactory.hpp"
#include "m2n/PointToPointComFactory.hpp"
#include "com/CommunicationFactory.hpp"

using namespace precice;
using namespace partition;

BOOST_AUTO_TEST_SUITE(PartitionTests)
BOOST_AUTO_TEST_SUITE(ProvidedBoundingBoxTests)

void setupParallelEnvironment(m2n::PtrM2N m2n){
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
    masterSlaveCom->requestConnection( "SolidMaster", "SolidSlaves", 0 , 2 );
  }
  else if(utils::Parallel::getProcessRank() == 3){//Slave2
    masterSlaveCom->requestConnection( "SolidMaster", "SolidSlaves", 1 , 2 );
  }
}


void setupM2NEnvironment(m2n::PtrM2N m2n){
  assertion(utils::Parallel::getCommunicatorSize() == 4);

  com::PtrCommunication masterSlaveCom =
        com::PtrCommunication(new com::MPIDirectCommunication());
  utils::MasterSlave::_communication = masterSlaveCom;

  utils::Parallel::synchronizeProcesses();

  if (utils::Parallel::getProcessRank() == 0){ //Master Fluid
    utils::Parallel::splitCommunicator( "FluidMaster" );
    m2n->acceptMasterConnection ( "FluidMaster", "SolidMaster");
    utils::MasterSlave::_rank = 0;
    utils::MasterSlave::_size = 2;
    utils::MasterSlave::_slaveMode = false;
    utils::MasterSlave::_masterMode = true;
   
  }
  else if(utils::Parallel::getProcessRank() == 1){//Slave1
    utils::Parallel::splitCommunicator( "FluidSlave");
    utils::MasterSlave::_rank = 1;
    utils::MasterSlave::_size = 2;
    utils::MasterSlave::_slaveMode = true;
    utils::MasterSlave::_masterMode = false;
   
  }
  else if(utils::Parallel::getProcessRank() == 2){//Master Solid
    utils::Parallel::splitCommunicator( "SolidMaster" );
    m2n->requestMasterConnection ( "FluidMaster", "SolidMaster" );
    utils::MasterSlave::_rank = 0;
    utils::MasterSlave::_size = 2;
    utils::MasterSlave::_slaveMode = false;
    utils::MasterSlave::_masterMode = true;
   
  }
  else if(utils::Parallel::getProcessRank() == 3){//Slave2
    utils::Parallel::splitCommunicator( "SolidSlave");
    utils::MasterSlave::_rank = 1;
    utils::MasterSlave::_size = 2;
    utils::MasterSlave::_slaveMode = true;
    utils::MasterSlave::_masterMode = false;
   
  }

  if(utils::Parallel::getProcessRank() == 0){//Master
    masterSlaveCom->acceptConnection ( "FluidMaster", "FluidSlave");
    masterSlaveCom->setRankOffset(1);
  }
  else if(utils::Parallel::getProcessRank() == 1){//Slave
    masterSlaveCom->requestConnection( "FluidMaster", "FluidSlave", 1 , 1 );
  }

  if(utils::Parallel::getProcessRank() == 2){//Master
    masterSlaveCom->acceptConnection ( "SolidMaster", "SolidSlave");
    masterSlaveCom->setRankOffset(1);
  }
  else if(utils::Parallel::getProcessRank() == 3){//Slave
    masterSlaveCom->requestConnection( "SolidMaster", "SolidSlave", 1 , 1 );
  }
}


void setupP2PEnvironment(m2n::PtrM2N m2n){assertion(utils::Parallel::getCommunicatorSize() == 4);

  com::PtrCommunication masterSlaveCom = com::PtrCommunication(new com::SocketCommunication());
  utils::MasterSlave::_communication = masterSlaveCom;

  utils::Parallel::synchronizeProcesses();

  if (utils::Parallel::getProcessRank() == 0){ //Master Fluid
    utils::Parallel::splitCommunicator( "FluidMaster" );
    m2n->acceptMasterConnection ( "FluidMaster", "SolidMaster");
    utils::MasterSlave::_rank = 0;
    utils::MasterSlave::_size = 2;
    utils::MasterSlave::_slaveMode = false;
    utils::MasterSlave::_masterMode = true;
   
  }
  else if(utils::Parallel::getProcessRank() == 1){//Slave1
    utils::Parallel::splitCommunicator( "FluidSlave");
    utils::MasterSlave::_rank = 1;
    utils::MasterSlave::_size = 2;
    utils::MasterSlave::_slaveMode = true;
    utils::MasterSlave::_masterMode = false;
   
  }
  else if(utils::Parallel::getProcessRank() == 2){//Master Solid
    utils::Parallel::splitCommunicator( "SolidMaster" );
    m2n->requestMasterConnection ( "FluidMaster", "SolidMaster" );
    utils::MasterSlave::_rank = 0;
    utils::MasterSlave::_size = 2;
    utils::MasterSlave::_slaveMode = false;
    utils::MasterSlave::_masterMode = true;
   
  }
  else if(utils::Parallel::getProcessRank() == 3){//Slave2
    utils::Parallel::splitCommunicator( "SolidSlave");
    utils::MasterSlave::_rank = 1;
    utils::MasterSlave::_size = 2;
    utils::MasterSlave::_slaveMode = true;
    utils::MasterSlave::_masterMode = false;
   
  }

  if(utils::Parallel::getProcessRank() == 0){//Master
    masterSlaveCom->acceptConnection ( "FluidMaster", "FluidSlave");
    masterSlaveCom->setRankOffset(1);
  }
  else if(utils::Parallel::getProcessRank() == 1){//Slave
    masterSlaveCom->requestConnection( "FluidMaster", "FluidSlave", 1 , 1 );
  }

  if(utils::Parallel::getProcessRank() == 2){//Master
    masterSlaveCom->acceptConnection ( "SolidMaster", "SolidSlave");
    masterSlaveCom->setRankOffset(1);
  }
  else if(utils::Parallel::getProcessRank() == 3){//Slave
    masterSlaveCom->requestConnection( "SolidMaster", "SolidSlave", 1 , 1 );
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


void process(std::vector<double> &data)
{
  for (auto &elem : data) {
    elem += utils::MasterSlave::_rank + 1;
  }
}


BOOST_AUTO_TEST_CASE(TestProvidedBoundingBox2D, * testing::OnSize(4))
{
  com::PtrCommunication participantCom =
      com::PtrCommunication(new com::MPIDirectCommunication());
  m2n::DistributedComFactory::SharedPointer distrFactory = m2n::DistributedComFactory::SharedPointer(
      new m2n::GatherScatterComFactory(participantCom));
  m2n::PtrM2N m2n = m2n::PtrM2N(new m2n::M2N(participantCom, distrFactory));

  setupParallelEnvironment(m2n);

  int dimensions = 2;
  bool flipNormals = true;

  BOOST_TEST(utils::Parallel::getCommunicatorSize()==4); 
  

  if (utils::Parallel::getProcessRank() != 0) { //NASTIN
  
    mesh::PtrMesh pSolidzMesh(new mesh::Mesh("SolidzMesh", dimensions, flipNormals));

    if(utils::Parallel::getProcessRank() == 1){//Master
      Eigen::VectorXd position(dimensions);
      position << -1.0, 0.0;
      pSolidzMesh->createVertex(position);
      position << 1.0, 2.0;
      pSolidzMesh->createVertex(position);
      position << 5.0, 3.0;
      pSolidzMesh->createVertex(position);
    }   
  
    else if(utils::Parallel::getProcessRank() == 2){//Slave1
      Eigen::VectorXd position(dimensions);
      position << 1.0, 3.5;
      mesh::Vertex& v3 = pSolidzMesh->createVertex(position);
      position << 0.0, 4.5;
      mesh::Vertex& v4 = pSolidzMesh->createVertex(position);
      pSolidzMesh->createEdge(v3,v4);
    }
    else if(utils::Parallel::getProcessRank() == 3){//Slave2
      Eigen::VectorXd position(dimensions);
      position << 2.5, 5.5;
      mesh::Vertex& v5 = pSolidzMesh->createVertex(position);
      position << 4.5, 7.0;
      mesh::Vertex& v6 = pSolidzMesh->createVertex(position); 
      pSolidzMesh->createEdge(v5,v6);
    }
 
    double safetyFactor = 0.0;
    bool hasToSend = true;
    mesh::Mesh::BoundingBoxMap send_globalBB;
    mesh::Mesh::BoundingBox localBB;
    pSolidzMesh->computeState();

    for (int i=0; i < 3; i++) {
      for (int j=0; j < dimensions; j++) {
        localBB.push_back(std::make_pair(-1,-1));
      }
      send_globalBB[i]=localBB;
      localBB.clear();
    }

    ProvidedBoundingBox part(pSolidzMesh, hasToSend, safetyFactor);
    part.setm2n(m2n);
    part.communicateBoundingBox();


    if(utils::Parallel::getProcessRank() == 1){//Master
      BOOST_TEST(part._globalBB.size()==3);
      BOOST_TEST(part._bb.size()==2);      
    } else if(utils::Parallel::getProcessRank() == 2)
    {
      BOOST_TEST(part._bb[1].first==3.5);
      BOOST_TEST(part._bb[1].second==4.5);
    }
    
  }
  else{//SOLIDZ

    
    mesh::PtrMesh pSolidzMesh(new mesh::Mesh("SolidzMesh", dimensions, flipNormals));    
    double safetyFactor = 0;
    Eigen::VectorXd position(dimensions);
      position << 0.0, 0.0;
      pSolidzMesh->createVertex(position);
      position << 0.0, 2.0;
      pSolidzMesh->createVertex(position);
      position << 0.0, 3.0;
      pSolidzMesh->createVertex(position);
   
    mesh::Mesh::BoundingBoxMap received_globalBB;    
    mesh::Mesh::BoundingBox localBB;
    pSolidzMesh->computeState();

    for (int i=0; i < 3; i++) {
      for (int j=0; j < dimensions; j++) {
        localBB.push_back(std::make_pair(-1,-1));
      }
      received_globalBB[i]=localBB;
      localBB.clear();
    }
    

    ReceivedBoundingBox part(pSolidzMesh, safetyFactor, ReceivedBoundingBox::FILTER_FIRST);
    part.setm2n(m2n);
    part.communicateBoundingBox();

    BOOST_TEST(part._bb[1].first==0);
    BOOST_TEST(part._bb[1].second==3);

    BOOST_TEST(part.remoteParComSize==3);    
    BOOST_TEST(part._globalBB.size()==3);

    BOOST_TEST(part._globalBB[0][0].first==-1.0);
    BOOST_TEST(part._globalBB[0][0].second==5.0);
    BOOST_TEST(part._globalBB[0][1].first==0.0);
    BOOST_TEST(part._globalBB[0][1].second==3.0);
    BOOST_TEST(part._globalBB[1][0].first==0.0);
    BOOST_TEST(part._globalBB[1][0].second==1.0);
    BOOST_TEST(part._globalBB[1][1].first==3.5);
    BOOST_TEST(part._globalBB[1][1].second==4.5);
    BOOST_TEST(part._globalBB[2][0].first==2.5);
    BOOST_TEST(part._globalBB[2][0].second==4.5);
    BOOST_TEST(part._globalBB[2][1].first==5.5);
    BOOST_TEST(part._globalBB[2][1].second==7.0);
    }

  tearDownParallelEnvironment();
  
}

BOOST_AUTO_TEST_CASE(TestProvidedBoundingBox3D, * testing::OnSize(4))
{
  com::PtrCommunication participantCom =
      com::PtrCommunication(new com::MPIDirectCommunication());
  m2n::DistributedComFactory::SharedPointer distrFactory = m2n::DistributedComFactory::SharedPointer(
      new m2n::GatherScatterComFactory(participantCom));
  m2n::PtrM2N m2n = m2n::PtrM2N(new m2n::M2N(participantCom, distrFactory));

  setupParallelEnvironment(m2n);

  int dimensions = 3;
  bool flipNormals = true;

  BOOST_TEST(utils::Parallel::getCommunicatorSize()==4); 
  

  if (utils::Parallel::getProcessRank() != 0) { //NASTIN
  
    mesh::PtrMesh pSolidzMesh(new mesh::Mesh("SolidzMesh", dimensions, flipNormals));

    if(utils::Parallel::getProcessRank() == 1){//Master
      Eigen::VectorXd position(dimensions);
      position << -1.0, 0.0, 1.0;
      pSolidzMesh->createVertex(position);
      position << 1.0, 2.0, -1.0;
      pSolidzMesh->createVertex(position);
      position << 5.0, 3.0, 0.0;
      pSolidzMesh->createVertex(position);
    }   
  
    else if(utils::Parallel::getProcessRank() == 2){//Slave1
      Eigen::VectorXd position(dimensions);
      position << 1.0, 3.5, 0.0;
      mesh::Vertex& v3 = pSolidzMesh->createVertex(position);
      position << 0.0, 4.5, -2.0;
      mesh::Vertex& v4 = pSolidzMesh->createVertex(position);
      pSolidzMesh->createEdge(v3,v4);
    }
    else if(utils::Parallel::getProcessRank() == 3){//Slave2
      Eigen::VectorXd position(dimensions);
      position << 2.5, 5.5, 1.5;
      mesh::Vertex& v5 = pSolidzMesh->createVertex(position);
      position << 4.5, 7.0, 3.5;
      mesh::Vertex& v6 = pSolidzMesh->createVertex(position); 
      pSolidzMesh->createEdge(v5,v6);
    }
 
    double safetyFactor = 0.0;
    bool hasToSend = true;
    mesh::Mesh::BoundingBoxMap send_globalBB;
    mesh::Mesh::BoundingBox localBB;
    pSolidzMesh->computeState();

    for (int i=0; i < 3; i++) {
      for (int j=0; j < dimensions; j++) {
        localBB.push_back(std::make_pair(-1,-1));
      }
      send_globalBB[i]=localBB;
      localBB.clear();
    }

    ProvidedBoundingBox part(pSolidzMesh, hasToSend, safetyFactor);
    part.setm2n(m2n);
    part.communicateBoundingBox();


    if(utils::Parallel::getProcessRank() == 1){//Master
      BOOST_TEST(part._globalBB.size()==3);
      BOOST_TEST(part._bb.size()==3);
    } 
    
  }
  else{//SOLIDZ

    
    mesh::PtrMesh pSolidzMesh(new mesh::Mesh("SolidzMesh", dimensions, flipNormals));    
    double safetyFactor = 0;
    Eigen::VectorXd position(dimensions);
    position << 0.0, 0.0, 0.0;
    pSolidzMesh->createVertex(position);
    position << 0.0, 2.0, -1.0;
    pSolidzMesh->createVertex(position);
    position << 0.0, 3.0, 5.0;
    pSolidzMesh->createVertex(position);
   
    mesh::Mesh::BoundingBoxMap received_globalBB;    
    mesh::Mesh::BoundingBox localBB;
    pSolidzMesh->computeState();

    for (int i=0; i < 3; i++) {
      for (int j=0; j < dimensions; j++) {
        localBB.push_back(std::make_pair(-1,-1));
      }
      received_globalBB[i]=localBB;
      localBB.clear();
    }
    

    ReceivedBoundingBox part(pSolidzMesh, safetyFactor, ReceivedBoundingBox::FILTER_FIRST);
    part.setm2n(m2n);
    part.communicateBoundingBox();

    BOOST_TEST(part._bb[1].first==0);
    BOOST_TEST(part._bb[1].second==3);

    BOOST_TEST(part.remoteParComSize==3);    
    BOOST_TEST(part._globalBB.size()==3);

    BOOST_TEST(part._globalBB[0][0].first==-1.0);
    BOOST_TEST(part._globalBB[0][0].second==5.0);
    BOOST_TEST(part._globalBB[0][1].first==0.0);
    BOOST_TEST(part._globalBB[0][1].second==3.0);
    BOOST_TEST(part._globalBB[0][2].first==-1.0);
    BOOST_TEST(part._globalBB[0][2].second==1.0);
    BOOST_TEST(part._globalBB[1][0].first==0.0);
    BOOST_TEST(part._globalBB[1][0].second==1.0);
    BOOST_TEST(part._globalBB[1][1].first==3.5);
    BOOST_TEST(part._globalBB[1][1].second==4.5);
    BOOST_TEST(part._globalBB[1][2].first==-2.0);
    BOOST_TEST(part._globalBB[1][2].second==0.0);
    BOOST_TEST(part._globalBB[2][0].first==2.5);
    BOOST_TEST(part._globalBB[2][0].second==4.5);
    BOOST_TEST(part._globalBB[2][1].first==5.5);
    BOOST_TEST(part._globalBB[2][1].second==7.0);
    BOOST_TEST(part._globalBB[2][2].first==1.5);
    BOOST_TEST(part._globalBB[2][2].second==3.5);

    }

  tearDownParallelEnvironment();
  
}

BOOST_AUTO_TEST_CASE(TestInitialCommunicationMap, * testing::OnSize(4))
{
  com::PtrCommunication participantCom =
      com::PtrCommunication(new com::MPIDirectCommunication());
  m2n::DistributedComFactory::SharedPointer distrFactory = m2n::DistributedComFactory::SharedPointer(
      new m2n::GatherScatterComFactory(participantCom));
  m2n::PtrM2N m2n = m2n::PtrM2N(new m2n::M2N(participantCom, distrFactory));

  setupParallelEnvironment(m2n);

  int dimensions = 3;
  bool flipNormals = true;

  BOOST_TEST(utils::Parallel::getCommunicatorSize()==4); 
  

  if (utils::Parallel::getProcessRank() != 0) { //NASTIN
  
    mesh::PtrMesh pSolidzMesh(new mesh::Mesh("SolidzMesh", dimensions, flipNormals));

    if(utils::Parallel::getProcessRank() == 1){//Master
      Eigen::VectorXd position(dimensions);
      position << -1.0, 0.0, 1.0;
      pSolidzMesh->createVertex(position);
      position << 1.0, 2.0, -1.0;
      pSolidzMesh->createVertex(position);
      position << 5.0, 3.0, 0.0;
      pSolidzMesh->createVertex(position);
    }   
  
    else if(utils::Parallel::getProcessRank() == 2){//Slave1
      Eigen::VectorXd position(dimensions);
      position << 1.0, 3.5, 0.0;
      mesh::Vertex& v3 = pSolidzMesh->createVertex(position);
      position << 0.0, 4.5, -2.0;
      mesh::Vertex& v4 = pSolidzMesh->createVertex(position);
      pSolidzMesh->createEdge(v3,v4);
    }
    else if(utils::Parallel::getProcessRank() == 3){//Slave2
      Eigen::VectorXd position(dimensions);
      position << 2.5, 5.5, 1.5;
      mesh::Vertex& v5 = pSolidzMesh->createVertex(position);
      position << 4.5, 7.0, 3.5;
      mesh::Vertex& v6 = pSolidzMesh->createVertex(position); 
      pSolidzMesh->createEdge(v5,v6);
    }
 
    double safetyFactor = 0.0;
    bool hasToSend = true;
    mesh::Mesh::BoundingBoxMap send_globalBB;
    mesh::Mesh::BoundingBox localBB;
    pSolidzMesh->computeState();

    for (int i=0; i < 3; i++) {
      for (int j=0; j < dimensions; j++) {
        localBB.push_back(std::make_pair(-1,-1));
      }
      send_globalBB[i]=localBB;
      localBB.clear();
    }

    ProvidedBoundingBox part(pSolidzMesh, hasToSend, safetyFactor);
    part.setm2n(m2n);
    part.computeBoundingBox();
    
    
    //test whether we receive correct feedbackmap via masters
    if(utils::Parallel::getProcessRank() == 1) {
      BOOST_TEST(part.remoteParComSize==3);
      BOOST_TEST(part.received_feedbackMap.size()==3);
      mesh::Mesh::FeedbackMap providedFeedbackMap;
      std::vector<int> test;
      for (int i=0; i < 3; i++) {
        for (int j=0; j < 3; j++) {
          if (i!=j) {
            test.push_back(j+1);
          }
        }
        providedFeedbackMap[i]=test;
        test.clear();
        BOOST_TEST(part.received_feedbackMap[i]==providedFeedbackMap[i]);
      }
    }

    // test whether initial communication map is built correctly
    if(utils::Parallel::getProcessRank() == 1) {
      BOOST_TEST(pSolidzMesh->getCommunicationMap().size()==2);
      BOOST_TEST(pSolidzMesh->getCommunicationMap()[1][0]==-1);
      BOOST_TEST(pSolidzMesh->getCommunicationMap()[2][0]==-1);
    } else if(utils::Parallel::getProcessRank() == 2) {
      BOOST_TEST(part.received_feedbackMap.size()==3);
      BOOST_TEST(pSolidzMesh->getCommunicationMap().size()==2);
      BOOST_TEST(pSolidzMesh->getCommunicationMap()[0][0]==-1);
      BOOST_TEST(pSolidzMesh->getCommunicationMap()[2][0]==-1);
    }else if(utils::Parallel::getProcessRank() == 3) {
      BOOST_TEST(part.received_feedbackMap.size()==3);
      BOOST_TEST(pSolidzMesh->getCommunicationMap().size()==2);
      BOOST_TEST(pSolidzMesh->getCommunicationMap()[0][0]==-1);
      BOOST_TEST(pSolidzMesh->getCommunicationMap()[1][0]==-1);
    }

    // test whether each rank has correct copy of vertexcounters 
    BOOST_TEST(part.vertexCounters[0]==3);
    BOOST_TEST(part.vertexCounters[1]==5);
    BOOST_TEST(part.vertexCounters[2]==7);
  }
  else{//SOLIDZ

    mesh::Mesh::FeedbackMap receivedFeedbackMap;
    std::vector<int> test;
    for (int i=0; i < 3; i++) {
      for (int j=0; j < 3; j++) {        
        if (i!=j) {
          test.push_back(j+1);                 
        }
      }
      receivedFeedbackMap[i]=test;
      test.clear();
    }
    m2n->getMasterCommunication()->send(3, 0);
    com::CommunicateBoundingBox(m2n->getMasterCommunication()).sendFeedbackMap(receivedFeedbackMap, 0);

  }

  tearDownParallelEnvironment();
  
}

BOOST_AUTO_TEST_CASE(P2PComMeshTest, * testing::OnSize(4))
{
  
  assertion(utils::Parallel::getCommunicatorSize() == 4);

//  com::PtrCommunicationFactory cf(new com::SocketCommunicationFactory);
  com::PtrCommunicationFactory cf(new com::MPIPortsCommunicationFactory);
  
  utils::MasterSlave::_communication = std::make_shared<com::MPIDirectCommunication>();

  mesh::PtrMesh mesh(new mesh::Mesh("Mesh", 2, true));

  int dimensions = 2;
  bool flipNormals = true;


  switch (utils::Parallel::getProcessRank()) {
  case 0: {
    utils::Parallel::splitCommunicator("Fluid.Master");

    utils::MasterSlave::_rank       = 0;
    utils::MasterSlave::_size       = 2;
    utils::MasterSlave::_masterMode = true;
    utils::MasterSlave::_slaveMode  = false;

    utils::MasterSlave::_communication->acceptConnection("Fluid.Master", "Fluid.Slave");
    utils::MasterSlave::_communication->setRankOffset(1);

    Eigen::VectorXd position(dimensions);
    position <<0.5, 0.0;
    mesh::Vertex& v1 = mesh->createVertex(position);
    position << 1.0, 2.0;
    mesh::Vertex& v2 = mesh->createVertex(position);
    mesh->createEdge(v1, v2);

    mesh->getCommunicationMap()[0].push_back(-1);

    
    break;
  }
  case 1: {
    utils::Parallel::splitCommunicator("Fluid.Slave");

    utils::MasterSlave::_rank       = 1;
    utils::MasterSlave::_size       = 2;
    utils::MasterSlave::_masterMode = false;
    utils::MasterSlave::_slaveMode  = true;

    utils::MasterSlave::_communication->requestConnection("Fluid.Master", "Fluid.Slave", 0, 1);

    Eigen::VectorXd position(dimensions);
    position <<0.5, 0.0;
    mesh::Vertex& v1 = mesh->createVertex(position);
    position << 1.0, 2.0;
    mesh::Vertex& v2 = mesh->createVertex(position);
    mesh->createEdge(v1, v2);

    mesh->getCommunicationMap()[0].push_back(-1);
    mesh->getCommunicationMap()[1].push_back(-1);



    break;
  }
  case 2: {
    utils::Parallel::splitCommunicator("Solid.Master");

    utils::MasterSlave::_rank       = 0;
    utils::MasterSlave::_size       = 2;
    utils::MasterSlave::_masterMode = true;
    utils::MasterSlave::_slaveMode  = false;

    utils::MasterSlave::_communication->acceptConnection("Solid.Master", "Solid.Slave");
    utils::MasterSlave::_communication->setRankOffset(1);


    Eigen::VectorXd position(dimensions);
    position << 1.5, 0.0;
    mesh::Vertex& v1 = mesh->createVertex(position);
    position << 1.0, 2.0;
    mesh::Vertex& v2 = mesh->createVertex(position);
    mesh->createEdge(v1, v2);

    mesh->getCommunicationMap()[0].push_back(-1);
    mesh->getCommunicationMap()[1].push_back(-1);
    

    break;
  }
  case 3: {
    utils::Parallel::splitCommunicator("Solid.Slave");

    utils::MasterSlave::_rank       = 1;
    utils::MasterSlave::_size       = 2;
    utils::MasterSlave::_masterMode = false;
    utils::MasterSlave::_slaveMode  = true;

    utils::MasterSlave::_communication->requestConnection("Solid.Master", "Solid.Slave", 0, 1);

    Eigen::VectorXd position(dimensions);
    position <<1.0, 0.0;
    mesh::Vertex& v1 = mesh->createVertex(position);
    position << 1.0, 2.0;
    mesh::Vertex& v2 = mesh->createVertex(position);
    mesh->createEdge(v1, v2);

    mesh->getCommunicationMap()[1].push_back(-1);


    break;
  }
  }

  m2n::PointToPointCommunication c(cf, mesh);

  if (utils::Parallel::getProcessRank() < 2) {
  
    c.requestConnection("Solid", "Fluid");
    c.sendMesh(*mesh);
    BOOST_TEST(mesh->getID()==0);
   
  } else {

    c.acceptConnection("Solid", "Fluid");
    c.receiveMesh(*mesh);
    BOOST_TEST(mesh->getID()==0);     

      if(utils::Parallel::getProcessRank() ==2 )
      {
          BOOST_TEST(mesh->vertices()[2].getCoords()[0]==0.50);
          BOOST_TEST(mesh->vertices().size()==6);
      }

      if(utils::Parallel::getProcessRank() ==3 )
      {
          BOOST_TEST(mesh->vertices()[2].getCoords()[0]==0.50);
          BOOST_TEST(mesh->vertices().size()==4);
      }


    
  }

  utils:: MasterSlave::_communication.reset();
  utils::MasterSlave::reset();

  utils::Parallel::synchronizeProcesses();
  utils::Parallel::clearGroups();
}




BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()
#endif // PRECICE_NO_MPI
