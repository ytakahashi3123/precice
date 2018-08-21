#pragma once

#include "Partition.hpp"
#include "logging/Logger.hpp"
#include <vector>
#include "mesh/Vertex.hpp"
#include "mesh/Mesh.hpp"


// Forward delcration to friend the boost test struct

namespace PartitionTests {
namespace ReceivedBoundingBoxTests {
struct TestReceivedBoundingBox2D;
struct TestReceivedBoundingBox3D;
}
}

namespace PartitionTests {
namespace ProvidedBoundingBoxTests {
struct TestProvidedBoundingBox;
}
}

namespace precice {
namespace partition {
/**
 * @brief A partition that is computed from a mesh received from another participant.
 *
 * A mesh is received by the master rank and re-partitioned among all slave ranks.
 * Afterwards necessary distribution data structures are set up.
 */
class ReceivedBoundingBox : public Partition
{
public:

   /// Constructor
  ReceivedBoundingBox (mesh::PtrMesh mesh, double safetyFactor, mesh::Mesh::BoundingBoxMap globalBB);
  virtual ~ReceivedBoundingBox() {}
  /// The mesh is received from another participant.
  virtual void communicate ();
  /// The mesh is re-partitioned and all distribution data structures are set up.
  virtual void compute ();
  void prepareBoundingBox();
  friend struct PartitionTests::ReceivedBoundingBoxTests::TestReceivedBoundingBox2D;
  friend struct PartitionTests::ReceivedBoundingBoxTests::TestReceivedBoundingBox3D;
  friend struct PartitionTests::ProvidedBoundingBoxTests::TestProvidedBoundingBox;
  bool CompareBoundingBox(mesh::Mesh::BoundingBox currentBB, mesh::Mesh::BoundingBox receivedBB); 

private:

  virtual void createOwnerInformation();
  void filterMesh(mesh::Mesh& filteredMesh, const bool filterByBB);
  mesh::Mesh::BoundingBoxMap _globalBB;
  mesh::Mesh::BoundingBox _bb;
  std::vector<int> feedback; 
  mesh::Mesh::FeedbackMap  feedbackMap; // int : each rank, vect: connected ranks
  int _dimensions;
  double _safetyFactor;
  int numberOfVertices;
  static logging::Logger _log;
};
}} // namespace precice, partition
