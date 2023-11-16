#pragma once
#include <wmtk/PolygonMesh.hpp>
#include <wmtk/operations/TupleOperation.hpp>
#include "PolygonMeshOperation.hpp"

namespace wmtk::operations {
namespace polygon_mesh {
class JoinFace;
}

template <>
struct OperationSettings<polygon_mesh::JoinFace>
{
    // TODO Determine if operation settings are necessary
};

namespace polygon_mesh {

/**
 * @class Operation to remove an edge and merge the adjacent faces into one.
 */
class JoinFace : public PolygonMeshOperation
{
public:
    // constructor for default factory pattern construction
    JoinFace(Mesh& m, const Tuple& h, const OperationSettings<JoinFace>& settings);

    std::string name() const override;

    using PolygonMeshOperation::hash_accessor;

    /**
     * @brief Check the precondition that the two faces of the edge are distinct.
     *
     * @return true if the precondition is satisfied
     * @return false otherwise
     */
    bool precondition();

protected:
    bool execute() override;

private:
    const OperationSettings<JoinFace>& m_settings;

    Tuple m_tuple;
};

} // namespace polygon_mesh
} // namespace wmtk::operations
