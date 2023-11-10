#pragma once
#include <wmtk/PolygonMesh.hpp>
#include <wmtk/operations/TupleOperation.hpp>
#include "AtomicOperation.hpp"
#include "PolygonMeshOperation.hpp"

namespace wmtk::operations {
namespace polygon_mesh {
class Splice;
}

template <>
struct OperationSettings<polygon_mesh::Splice>
{
    // TODO Determine if operation settings are necessary
};

namespace polygon_mesh {

/**
 * @class Atomic operation for changing connectivity without changing the number of edges.
 */
class Splice : public AtomicOperation
{
public:
    // constructor for default factory pattern construction
    Splice(Mesh& m, const Tuple& h, const Tuple& g, const OperationSettings<Splice>& settings);

    std::string name() const override;

    /**
     * @brief Check the precondition that the mesh remains manifold after the splice.
     *
     * Splice always produces a manifold connectivity if there are no hole faces, but nonmanifold
     * vertices are possible on a mesh with boundary.
     *
     * @return true if the mesh remains manifold
     * @return false otherwise
     */
    bool precondition() const;

    /**
     * @brief Determine if the splice preserves the topology of the original mesh
     *
     * @return true if the splice will preserve the topology
     * @return false otherwise
     */
    bool is_topology_preserving() const;

protected:
    bool execute() override;

private:
    Tuple m_first_tuple;
    Tuple m_second_tuple;
    const OperationSettings<Splice>& m_settings;

    /**
     * @brief Helper function for determining the next halfedge after the splice is complete.
     *
     * @param h: halfedge index
     * @return next halfedge of h after the splice is complete
     */
    long next_after_splice(long h) const;
};

} // namespace polygon_mesh
} // namespace wmtk::operations
