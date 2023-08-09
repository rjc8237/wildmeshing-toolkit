#pragma once
#include <optional>
#include "Operation.hpp"

namespace wmtk {
class TriMeshSplitEdgeOperation : public Operation
{
public:
    TriMeshSplitEdgeOperation(
        Mesh& m,
        const Tuple& t,
        const OperationSettings<TriMeshSplitEdgeOperation> = {});

    std::string name() const override;


    std::vector<Tuple> triangle_onering() const;
    std::vector<Tuple> triangle_tworing() const;
    std::vector<Tuple> edge_onering() const;

    Tuple new_vertex() const;
    Tuple return_tuple() const;
    // std::vector<Tuple> new_triangles() const ;
    // std::vector<Tuple> new_edges() const ;

    // std::array<Tuple,2> spline_edges() const;
    // std::vector<Tuple> rib_edges() const;

    static PrimitiveType primitive_type() { return PrimitiveType::Edge; }

protected:
    bool execute() override;
    bool before() const override;

private:
    Tuple m_input_tuple;
    Tuple m_output_tuple;
};


} // namespace wmtk
