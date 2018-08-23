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

  /// Defines the typ of geometric filter used
  enum GeometricFilter {
    // @brief undefined
    UNDEFINED,
    // @brief No geometric filter used (e.g. for RBF mappings)
    NO_FILTER,
    // @brief Filter at master and communicate only filtered mesh.
    FILTER_FIRST,
    // @brief Broadcast first and filter then
    BROADCAST_FILTER
  };

   /// Constructor
  ReceivedBoundingBox (mesh::PtrMesh mesh, double safetyFactor, GeometricFilter geometricFilter);
  virtual ~ReceivedBoundingBox() {}
  /// The mesh is received from another participant.
  virtual void communicate ();
  /// The mesh is re-partitioned and all distribution data structures are set up.
  virtual void compute ();
  virtual void communicatePartition();
  virtual void computePartition();  
  friend struct PartitionTests::ReceivedBoundingBoxTests::TestReceivedBoundingBox2D;
  friend struct PartitionTests::ReceivedBoundingBoxTests::TestReceivedBoundingBox3D;
  friend struct PartitionTests::ProvidedBoundingBoxTests::TestProvidedBoundingBox;

private:

  
  bool CompareBoundingBox(mesh::Mesh::BoundingBox currentBB, mesh::Mesh::BoundingBox receivedBB);

  void filterMesh(mesh::Mesh& filteredMesh, const bool filterByBB);

  void prepareBoundingBox();

  bool isVertexInBB(const mesh::Vertex& vertex);

  virtual void createOwnerInformation();

  /// Helper function for 'createOwnerFunction' to set local owner information
  void setOwnerInformation(const std::vector<int> &ownerVec);
 
  mesh::Mesh::BoundingBoxMap _globalBB;
  mesh::Mesh::BoundingBox _bb;
  std::vector<int> feedback; 
  mesh::Mesh::FeedbackMap  feedbackMap; // int : each rank, vect: connected ranks
  int _dimensions;
  double _safetyFactor;
  int numberOfVertices;
  static logging::Logger _log;
  std::vector<int> vertexCounters;
  GeometricFilter _geometricFilter;
  mesh::Mesh::FeedbackMap localCommunicationMap;
};
}} // namespace precice, partition
