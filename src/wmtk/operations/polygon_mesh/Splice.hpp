#pragma once
#include <wmtk/operations/TupleOperation.hpp>
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
class Splice : public PolygonMeshOperation, private TupleOperation
{
public:
    // constructor for default factory pattern construction
    Splice(Mesh& m, const Tuple& h, const Tuple& g, const OperationSettings<Splice>& settings);
    Splice(PolygonMesh& m, const Tuple& h, const Tuple& g, const OperationSettings<Splice>& settings);

    std::string name() const override;

    std::vector<Tuple> modified_primitives(PrimitiveType) const override;

    using PolygonMeshOperation::hash_accessor;

protected:
    bool execute() override;

private:
    const OperationSettings<Splice>& m_settings;
};

} // namespace polygon_mesh
} // namespace wmtk::operations
