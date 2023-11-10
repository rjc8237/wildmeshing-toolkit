#pragma once
#include <wmtk/PolygonMesh.hpp>
#include <wmtk/operations/TupleOperation.hpp>
#include "PolygonMeshOperation.hpp"

namespace wmtk::operations {
namespace polygon_mesh {
class SplitVertex;
}

template <>
struct OperationSettings<polygon_mesh::SplitVertex>
{
    // TODO Determine if operation settings are necessary
};

namespace polygon_mesh {

/**
 * @class Operation to split the vertex incident to h and g into two vertices and connect
 * them with a new edge.
 */
class SplitVertex : public PolygonMeshOperation
{
public:
    // constructor for default factory pattern construction
    SplitVertex(
        Mesh& m,
        const Tuple& h,
        const Tuple& g,
        const OperationSettings<SplitVertex>& settings);

    std::string name() const override;

    using PolygonMeshOperation::hash_accessor;

    /**
     * @brief Check the precondition that the halfedges are distinct and point to the same vertex.
     *
     * @return true if the precondition is satisfied
     * @return false otherwise
     */
    bool precondition();

    /**
     * @brief Return tuple corresponding to the halfedge created in the same split face as g.
     *
     * The vertex of the tuple is the base vertex of the halfedge.
     */
    Tuple return_tuple() const;

protected:
    bool execute() override;

private:
    const OperationSettings<SplitVertex>& m_settings;

    Tuple m_first_tuple;
    Tuple m_second_tuple;
    Tuple m_output_tuple;
};

} // namespace polygon_mesh
} // namespace wmtk::operations
