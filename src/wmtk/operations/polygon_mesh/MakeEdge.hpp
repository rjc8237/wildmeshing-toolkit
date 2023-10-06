#pragma once
#include <wmtk/operations/TupleOperation.hpp>
#include "PolygonMeshOperation.hpp"

namespace wmtk::operations {
namespace polygon_mesh {
class MakeEdge;
}

template <>
struct OperationSettings<polygon_mesh::MakeEdge>
{
    // TODO Determine if operation settings are necessary
};

namespace polygon_mesh {

/**
 * @class Atomic operation for adding a bubble component with a single edge.
 */
class MakeEdge : public PolygonMeshOperation, private TupleOperation
{
public:
    // constructor for default factory pattern construction
    MakeEdge(Mesh& m, const OperationSettings<MakeEdge>& settings);
    MakeEdge(PolygonMesh& m, const OperationSettings<MakeEdge>& settings);

    std::string name() const override;

    /**
     * @brief Return tuple corresponding to the created edge
     */
    Tuple return_tuple() const;

    std::vector<Tuple> modified_primitives(PrimitiveType) const override;

    using PolygonMeshOperation::hash_accessor;

protected:
    Tuple m_output_tuple;
    bool execute() override;

private:
    const OperationSettings<MakeEdge>& m_settings;
};

} // namespace polygon_mesh
} // namespace wmtk::operations
