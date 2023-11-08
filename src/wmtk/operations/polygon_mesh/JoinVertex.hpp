#pragma once
#include <wmtk/PolygonMesh.hpp>
#include <wmtk/operations/TupleOperation.hpp>
#include "PolygonMeshOperation.hpp"

namespace wmtk::operations {
namespace polygon_mesh {
class JoinVertex;
}

template <>
struct OperationSettings<polygon_mesh::JoinVertex>
{
    // TODO Determine if operation settings are necessary
};

namespace polygon_mesh {

/**
 * @class Operation to remove an edge and merge the tip and tail into one.
 *
 * If the tip and tail are the same, just removes the edge.
 */
class JoinVertex : public PolygonMeshOperation
{
public:
    // constructor for default factory pattern construction
    JoinVertex(Mesh& m, const Tuple& h, const OperationSettings<JoinVertex>& settings);

    std::string name() const override;

    using PolygonMeshOperation::hash_accessor;

    /**
     * @brief Check the precondition that one of the two vertices has valence at least 2
     * and that distinct vertices are not both boundaries unless they are on a boundary edge.
     *
     * @return true if the precondition is satisfied
     * @return false otherwise
     */
    bool precondition();

protected:
    bool execute() override;

private:
    const OperationSettings<JoinVertex>& m_settings;

    Tuple m_tuple;
};

} // namespace polygon_mesh
} // namespace wmtk::operations
