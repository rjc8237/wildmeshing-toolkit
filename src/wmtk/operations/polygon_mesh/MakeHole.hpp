#pragma once
#include <wmtk/operations/TupleOperation.hpp>
#include "AtomicOperation.hpp"
#include "PolygonMeshOperation.hpp"

namespace wmtk::operations {
namespace polygon_mesh {
class MakeHole;
}

template <>
struct OperationSettings<polygon_mesh::MakeHole>
{
};

namespace polygon_mesh {

/**
 * @class Atomic operation for making a face of a mesh a hole.
 */
class MakeHole : public AtomicOperation
{
public:
    // constructor for default factory pattern construction
    MakeHole(Mesh& m, const Tuple& f, const OperationSettings<MakeHole>& settings);

    std::string name() const override;

    using PolygonMeshOperation::hash_accessor;

    /**
     * @brief Check the precondition that the mesh remains manifold after the face is made a hole.
     *
     * Check that making the face a hole does not make a vertex nonmanifold
     *
     * @return true if the mesh remains manifold
     * @return false otherwise
     */
    bool precondition();

protected:
    bool execute() override;
    bool precondition_at_vertex(long halfedge_id);

private:
    const OperationSettings<MakeHole>& m_settings;
    Tuple m_tuple;
};

} // namespace polygon_mesh
} // namespace wmtk::operations
