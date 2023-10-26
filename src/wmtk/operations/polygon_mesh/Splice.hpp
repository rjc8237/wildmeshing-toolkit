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

    using PolygonMeshOperation::hash_accessor;

protected:
    bool execute() override;

private:
    const OperationSettings<Splice>& m_settings;

    Tuple m_first_tuple;
    Tuple m_second_tuple;

    bool is_manifold_after_splice();
    bool is_hole(long face_id) const;
    long next_after_splice(long h);
};

} // namespace polygon_mesh
} // namespace wmtk::operations
