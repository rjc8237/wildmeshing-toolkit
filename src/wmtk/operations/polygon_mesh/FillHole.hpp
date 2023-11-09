#pragma once
#include <wmtk/operations/TupleOperation.hpp>
#include "AtomicOperation.hpp"
#include "PolygonMeshOperation.hpp"

namespace wmtk::operations {
namespace polygon_mesh {
class FillHole;
}

template <>
struct OperationSettings<polygon_mesh::FillHole>
{
};

namespace polygon_mesh {

/**
 * @class Atomic operation for filling the face of a mesh.
 */
class FillHole : public AtomicOperation
{
public:
    // constructor for default factory pattern construction
    FillHole(Mesh& m, const Tuple& f, const OperationSettings<FillHole>& settings);

    std::string name() const override;

    /**
     * @brief Check the precondition that the mesh remains manifold after the face is filled.
     *
     * Check that filling the face does not make a vertex nonmanifold
     *
     * @return true if the mesh remains manifold
     * @return false otherwise
     */
    bool precondition() const;

protected:
    bool execute() override;
    bool precondition_at_vertex(long halfedge_id) const;

private:
    const OperationSettings<FillHole>& m_settings;
    Tuple m_tuple;
};

} // namespace polygon_mesh
} // namespace wmtk::operations
