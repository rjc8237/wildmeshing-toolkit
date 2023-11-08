#pragma once
#include <wmtk/PolygonMesh.hpp>
#include <wmtk/operations/TupleOperation.hpp>
#include "PolygonMeshOperation.hpp"

namespace wmtk::operations {
namespace polygon_mesh {
class SplitFace;
}

template <>
struct OperationSettings<polygon_mesh::SplitFace>
{
    // TODO Determine if operation settings are necessary
};

namespace polygon_mesh {

/**
 * @class Operation to split the face incident to h and g into two faces with a new diagonal
 *  between the two vertices denoted by h and g respectively.
 *
 * The added face is always non-boundary; the face to be split may be boundary or not.
 * TODO: Check this statement.
 */
class SplitFace : public PolygonMeshOperation
{
public:
    // constructor for default factory pattern construction
    SplitFace(
        Mesh& m,
        const Tuple& h,
        const Tuple& g,
        const OperationSettings<SplitFace>& settings);

    std::string name() const override;

    using PolygonMeshOperation::hash_accessor;

    /**
     * @brief Check the precondition that the halfedges are distinct and belong to the same face.
     *
     * @return true if the precondition is satisfied
     * @return false otherwise
     */
    bool precondition();

protected:
    bool execute() override;

private:
    const OperationSettings<SplitFace>& m_settings;

    Tuple m_first_tuple;
    Tuple m_second_tuple;
};

} // namespace polygon_mesh
} // namespace wmtk::operations
