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
#include "mapping/NearestNeighborMapping.hpp"
#include "mapping/NearestProjectionMapping.hpp"
#include "mapping/PetRadialBasisFctMapping.hpp"
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
    utils::Parallel::splitCommunicator( "Fluid");
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
    utils::Parallel::splitCommunicator( "Solid");
    utils::MasterSlave::_rank = 1;
    utils::MasterSlave::_size = 2;
    utils::MasterSlave::_slaveMode = true;
    utils::MasterSlave::_masterMode = false;
   
 }

  if(utils::Parallel::getProcessRank() == 0){//Master
    masterSlaveCom->acceptConnection ( "FluidMaster", "Fluid");
    masterSlaveCom->setRankOffset(1);
  }
  else if(utils::Parallel::getProcessRank() == 1){//Slave
    masterSlaveCom->requestConnection( "FluidMaster", "Fluid", 0 , 1 );
  }

  if(utils::Parallel::getProcessRank() == 2){//Master
    masterSlaveCom->acceptConnection ( "SolidMaster", "Solid");
    masterSlaveCom->setRankOffset(1);
  }
  else if(utils::Parallel::getProcessRank() == 3){//Slave
    masterSlaveCom->requestConnection( "SolidMaster", "Solid", 0 , 1 );
  }
}


void setupP2PEnvironment(m2n::PtrM2N m2n)
{
  assertion(utils::Parallel::getCommunicatorSize() == 4);

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

void setupMasterSlave( com::PtrCommunication masterSlaveCom){
  assertion(utils::Parallel::getCommunicatorSize() == 4);
        
  utils::MasterSlave::_communication = masterSlaveCom;

  utils::Parallel::synchronizeProcesses();

  if (utils::Parallel::getProcessRank() == 0){ //Master Fluid
    utils::Parallel::splitCommunicator( "FluidMaster" );
    utils::MasterSlave::_rank = 0;
    utils::MasterSlave::_size = 2;
    utils::MasterSlave::_slaveMode = false;
    utils::MasterSlave::_masterMode = true;
   
  }
  else if(utils::Parallel::getProcessRank() == 1){//Slave1
    utils::Parallel::splitCommunicator( "Fluid");
    utils::MasterSlave::_rank = 1;
    utils::MasterSlave::_size = 2;
    utils::MasterSlave::_slaveMode = true;
    utils::MasterSlave::_masterMode = false;
   
  }
  else if(utils::Parallel::getProcessRank() == 2){//Master Solid
    utils::Parallel::splitCommunicator( "SolidMaster" );
    utils::MasterSlave::_rank = 0;
    utils::MasterSlave::_size = 2;
    utils::MasterSlave::_slaveMode = false;
    utils::MasterSlave::_masterMode = true;
   
  }
  else if(utils::Parallel::getProcessRank() == 3){//Slave2
    utils::Parallel::splitCommunicator( "Solid");
    utils::MasterSlave::_rank = 1;
    utils::MasterSlave::_size = 2;
    utils::MasterSlave::_slaveMode = true;
    utils::MasterSlave::_masterMode = false;
   
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
    part.setM2N(m2n);
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
    mesh::PtrMesh pSolidzMesh_received(new mesh::Mesh("SolidzMesh", dimensions, flipNormals));
    mapping::PtrMapping boundingFromMapping = mapping::PtrMapping(new mapping::NearestNeighborMapping(mapping::Mapping::CONSISTENT, dimensions));
    mapping::PtrMapping boundingToMapping = mapping::PtrMapping(new mapping::NearestNeighborMapping(mapping::Mapping::CONSERVATIVE, dimensions));
    boundingFromMapping->setMeshes(pSolidzMesh_received, pSolidzMesh);
    boundingToMapping->setMeshes(pSolidzMesh, pSolidzMesh_received);
    
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
    

    ReceivedBoundingBox part(pSolidzMesh_received, safetyFactor, ReceivedBoundingBox::FILTER_FIRST);
    part.setM2N(m2n);
    part.setFromMapping(boundingFromMapping);
    part.setToMapping(boundingToMapping);
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
    part.setM2N(m2n);
    part.communicateBoundingBox();


    if(utils::Parallel::getProcessRank() == 1){//Master
      BOOST_TEST(part._globalBB.size()==3);
      BOOST_TEST(part._bb.size()==3);
    } 
    
  }
  else{//SOLIDZ

    
    mesh::PtrMesh pSolidzMesh(new mesh::Mesh("SolidzMesh", dimensions, flipNormals));
    mesh::PtrMesh pSolidzMesh_received(new mesh::Mesh("SolidzMesh", dimensions, flipNormals));
    mapping::PtrMapping boundingFromMapping = mapping::PtrMapping(new mapping::NearestNeighborMapping(mapping::Mapping::CONSISTENT, dimensions));
    mapping::PtrMapping boundingToMapping = mapping::PtrMapping(new mapping::NearestNeighborMapping(mapping::Mapping::CONSERVATIVE, dimensions));
    boundingFromMapping->setMeshes(pSolidzMesh_received, pSolidzMesh);
    boundingToMapping->setMeshes(pSolidzMesh, pSolidzMesh_received);

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
    

    ReceivedBoundingBox part(pSolidzMesh_received, safetyFactor, ReceivedBoundingBox::FILTER_FIRST);
    part.setM2N(m2n);
    part.setFromMapping(boundingFromMapping);
    part.setToMapping(boundingToMapping);

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
    part.setM2N(m2n);
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


BOOST_AUTO_TEST_CASE(TestM2NMeshExchange, * testing::OnSize(4))
{


//mesh creation
  int dimensions = 2;
  bool flipNormals = true;
  double safetyFactor = 0.1;
  bool hasToSend=true;
  mesh::PtrMesh mesh(new mesh::Mesh("mesh", dimensions, flipNormals));
  mesh::PtrMesh mesh_received(new mesh::Mesh("SolidzMesh", dimensions, flipNormals));  
  
  switch (utils::Parallel::getProcessRank()) {
  case 0: {
    Eigen::VectorXd position(dimensions);
    position <<0.5, 0.0;
    mesh::Vertex& v1 = mesh->createVertex(position);
    position << 1.5, 0.0;
    mesh::Vertex& v2 = mesh->createVertex(position);
    position <<2.0, 1.0;
    mesh::Vertex& v3 = mesh->createVertex(position);
    position << 0.5, 1.0;
    mesh::Vertex& v4 = mesh->createVertex(position);
    mesh->createEdge(v1, v2);
    mesh->createEdge(v2, v3);
    mesh->createEdge(v3, v4);
    mesh->createEdge(v4, v1);
    
    
    break;
  }
  case 1: {
    Eigen::VectorXd position(dimensions);
    position <<2.5, 0.0;
    mesh::Vertex& v1 = mesh->createVertex(position);
    position << 3.5, 0.0;
    mesh::Vertex& v2 = mesh->createVertex(position);
    position <<3.5, 1.0;
    mesh::Vertex& v3 = mesh->createVertex(position);
    position << 2.0, 1.0;
    mesh::Vertex& v4 = mesh->createVertex(position);
    mesh->createEdge(v1, v2);
    mesh->createEdge(v2, v3);
    mesh->createEdge(v3, v4);
    mesh->createEdge(v4, v1);

    break;
  }
  case 2: {
    Eigen::VectorXd position(dimensions);
    position <<0.5, 0.0;
    mesh::Vertex& v1 = mesh->createVertex(position);
    position << 2.0, 0.0;
    mesh::Vertex& v2 = mesh->createVertex(position);
    position <<2.0, -1.0;
    mesh::Vertex& v3 = mesh->createVertex(position);
    position << 0.5, -1.0;
    mesh::Vertex& v4 = mesh->createVertex(position);
    mesh->createEdge(v1, v2);
    mesh->createEdge(v2, v3);
    mesh->createEdge(v3, v4);
    mesh->createEdge(v4, v1);

    break;
  }
  case 3: {
    Eigen::VectorXd position(dimensions);
    position <<2.0, 0.0;
    mesh::Vertex& v1 = mesh->createVertex(position);
    position << 3.5, 0.0;
    mesh::Vertex& v2 = mesh->createVertex(position);
    position <<3.5, -1.0;
    mesh::Vertex& v3 = mesh->createVertex(position);
    position << 2.0, -1.0;
    mesh::Vertex& v4 = mesh->createVertex(position);
    mesh->createEdge(v1, v2);
    mesh->createEdge(v2, v3);
    mesh->createEdge(v3, v4);
    mesh->createEdge(v4, v1);

    break;
  }
  }

    mesh->computeState();
  
// create communicatror for master com and bb exchange/com/initial com_map

  com::PtrCommunication participantCom1 = com::PtrCommunication(new com::SocketCommunication());
  m2n::DistributedComFactory::SharedPointer distrFactory = m2n::DistributedComFactory::SharedPointer(new m2n::GatherScatterComFactory(participantCom1));
  m2n::PtrM2N m2n = m2n::PtrM2N(new m2n::M2N(participantCom1, distrFactory));
  setupM2NEnvironment(m2n);
  

  com::PtrCommunication participantsCom =  com::PtrCommunication(new com::SocketCommunication());
  com::PtrCommunicationFactory participantComFactory =  com::PtrCommunicationFactory(new com::SocketCommunicationFactory);
  m2n::DistributedComFactory::SharedPointer distributionFactory = m2n::DistributedComFactory::SharedPointer(new m2n::PointToPointComFactory(participantComFactory));
  m2n::PtrM2N p2p = m2n::PtrM2N(new m2n::M2N(participantsCom, distributionFactory));
  //setupMasterSlave(participantsCom);
  

  if(utils::Parallel::getProcessRank() < 2)
  {
    p2p->createDistributedCommunication(mesh);
    ProvidedBoundingBox part(mesh, hasToSend, safetyFactor);
    part.setM2N(m2n);
    
    part.communicateBoundingBox();
    part.computeBoundingBox();

    BOOST_TEST(part._mesh->vertices().size()==4);
    BOOST_TEST(part.received_feedbackMap[0][0]==0);
    BOOST_TEST(part.received_feedbackMap[0][1]==1);
    BOOST_TEST(part.received_feedbackMap[1][0]==0);
    BOOST_TEST(part.received_feedbackMap[1][1]==1);
    BOOST_TEST(part.vertexCounters[0]==4);
    BOOST_TEST(part.vertexCounters[1]==8);

    BOOST_TEST(part.received_feedbackMap.size()==2);
    BOOST_TEST(part.received_feedbackMap[0][0]==0);
    BOOST_TEST(part.received_feedbackMap[0][1]==1);
    BOOST_TEST(part.received_feedbackMap[1][0]==0);
    BOOST_TEST(part.received_feedbackMap[1][1]==1);

    BOOST_TEST(part._mesh->getCommunicationMap().size()==2);
    BOOST_TEST(part._mesh->getCommunicationMap()[0][0]==-1);
    BOOST_TEST(part._mesh->getCommunicationMap()[1][0]==-1);
    BOOST_TEST(mesh->getCommunicationMap().size()==2);
    BOOST_TEST(mesh->getCommunicationMap()[0][0]==-1);
    BOOST_TEST(mesh->getCommunicationMap()[1][0]==-1);

    if(utils::Parallel::getProcessRank() == 0 )
    {
      BOOST_TEST(part._mesh->vertices()[0].getGlobalIndex()==0);
      BOOST_TEST(part._mesh->vertices()[1].getGlobalIndex()==1);
      BOOST_TEST(part._mesh->vertices()[2].getGlobalIndex()==2);
      BOOST_TEST(part._mesh->vertices()[3].getGlobalIndex()==3);
    }
    else
    {
      BOOST_TEST(part._mesh->vertices()[0].getGlobalIndex()==4);
      BOOST_TEST(part._mesh->vertices()[1].getGlobalIndex()==5);
      BOOST_TEST(part._mesh->vertices()[2].getGlobalIndex()==6);
      BOOST_TEST(part._mesh->vertices()[3].getGlobalIndex()==7);
    }

    
    part.setM2N(p2p);
    p2p->requestSlavesPreConnection("Solid", "Fluid");

    part.communicate();
    part.compute();

    //If some one need the list of vertices in CommunicationMap:
    
    // if(utils::Parallel::getProcessRank() == 1)
    //   {
    //     std::cout << "for global rank  " <<  utils::Parallel::getProcessRank() << std::endl;     
    //     for (auto & ranks : part._mesh->getCommunicationMap() ) {
    //       std::cout << "vertices for rank = " <<  ranks.first << "are " ;
    //       for (auto & vertices : ranks.second ) {
    //         std::cout << vertices << ", ";
    //       }
    //       std::cout <<  std::endl;
    //     }
    //   } 
  }
  else
  {
    p2p->createDistributedCommunication(mesh_received);
    mapping::PtrMapping boundingFromMapping = mapping::PtrMapping(new mapping::NearestNeighborMapping(mapping::Mapping::CONSISTENT, dimensions));
    mapping::PtrMapping boundingToMapping = mapping::PtrMapping(new mapping::NearestNeighborMapping(mapping::Mapping::CONSERVATIVE, dimensions));
    boundingFromMapping->setMeshes(mesh_received, mesh);
    boundingToMapping->setMeshes(mesh, mesh_received);

    ReceivedBoundingBox part(mesh_received, safetyFactor, ReceivedBoundingBox::BROADCAST_FILTER);
    
    part.setM2N(m2n);
    part.setFromMapping(boundingFromMapping);
    part.setToMapping(boundingToMapping);
    part.communicateBoundingBox();    
    part.computeBoundingBox();
    mesh_received->getVertexDistribution()[0].push_back(1);
    mesh_received->getVertexDistribution()[1].push_back(1);                                                   
    mesh_received->getCommunicationMap()[0].push_back(-1);
    mesh_received->getCommunicationMap()[1].push_back(-1);
    part.setM2N(p2p);
    p2p->acceptSlavesPreConnection("Solid", "Fluid");

    part.communicate();
    for (auto & vertex : part._mesh->vertices() ) {
      if (vertex.getGlobalIndex()==0) {
        BOOST_TEST(vertex.getCoords()[0]== 0.5);
        BOOST_TEST(vertex.getCoords()[1]== 0.0);
      } else if (vertex.getGlobalIndex()==1) {
        BOOST_TEST(vertex.getCoords()[0]== 1.5);
        BOOST_TEST(vertex.getCoords()[1]== 0.0);
      } else if (vertex.getGlobalIndex()==2) {
        BOOST_TEST(vertex.getCoords()[0]== 2.0);
        BOOST_TEST(vertex.getCoords()[1]== 1.0);
      } else if (vertex.getGlobalIndex()==3) {
        BOOST_TEST(vertex.getCoords()[0]== 0.5);
        BOOST_TEST(vertex.getCoords()[1]== 1.0);
      } else if (vertex.getGlobalIndex()==4) {
        BOOST_TEST(vertex.getCoords()[0]== 2.5);
        BOOST_TEST(vertex.getCoords()[1]== 0.0);
      } else if (vertex.getGlobalIndex()==5) {
        BOOST_TEST(vertex.getCoords()[0]== 3.5);
        BOOST_TEST(vertex.getCoords()[1]== 0.0);
      } else if (vertex.getGlobalIndex()==6) {
        BOOST_TEST(vertex.getCoords()[0]== 3.5);
        BOOST_TEST(vertex.getCoords()[1]== 1.0);
      } else if (vertex.getGlobalIndex()==7) {
        BOOST_TEST(vertex.getCoords()[0]== 2.0);
        BOOST_TEST(vertex.getCoords()[1]== 1.0);
      }
    }

      BOOST_TEST(part._mesh->vertices().size()==8);
      
      part.compute();
      
      BOOST_TEST(part._mesh->vertices().size()==2);
      if(utils::Parallel::getProcessRank() == 2)
      {
        BOOST_TEST(part.localCommunicationMap[1].size() == 0);
        BOOST_TEST(part.localCommunicationMap[0][0] == 0);
        BOOST_TEST(part.localCommunicationMap[0][1] == 1);
      } else
      {
        BOOST_TEST(part.localCommunicationMap[0].size() == 0);
        BOOST_TEST(part.localCommunicationMap[1][0] == 4);
        BOOST_TEST(part.localCommunicationMap[1][1] == 5);
      }

      

      //If some one need the list of vertices in CommunicationMap:
      
      // if(utils::Parallel::getProcessRank() == 2)
      // {
      //   std::cout << "for global rank  " <<  utils::Parallel::getProcessRank() << std::endl;     
      //   for (auto & ranks : part.localCommunicationMap ) {
      //     std::cout << "vertices for rank = " <<  ranks.first << "are " ;
      //     for (auto & vertices : ranks.second ) {
      //       std::cout << vertices << ", ";
      //     }
      //     std::cout <<  std::endl;
      //   }
      // }      
  } 
  tearDownParallelEnvironment();
}
  





BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()
#endif // PRECICE_NO_MPI
