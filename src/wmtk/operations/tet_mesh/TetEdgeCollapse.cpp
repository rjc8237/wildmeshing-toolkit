#include "TetEdgeCollapse.hpp"
#include <spdlog/spdlog.h>
#include <wmtk/SimplicialComplex.hpp>
#include <wmtk/TetMesh.hpp>
#include <wmtk/invariants/ValidTupleInvariant.hpp>
#include <wmtk/invariants/find_invariant_in_collection_by_type.hpp>

namespace wmtk::operations {

OperationSettings<tet_mesh::TetEdgeCollapse>::OperationSettings() {}

void OperationSettings<tet_mesh::TetEdgeCollapse>::initialize_invariants(const TetMesh& m)
{
    // outdated + is valid tuple
    invariants = basic_invariant_collection(m);
}

bool OperationSettings<tet_mesh::TetEdgeCollapse>::are_invariants_initialized() const
{
    return true;
}

namespace tet_mesh {

TetEdgeCollapse::TetEdgeCollapse(
    TetMesh& m,
    const Tuple& t,
    const OperationSettings<TetEdgeCollapse>& settings)
    : TetMeshOperation(m)
    , TupleOperation(settings.invariants, t)
// , m_settings(settings)
{
    assert(settings.are_invariants_initialized());
}

bool TetEdgeCollapse::execute()
{
    auto return_data = mesh().collapse_edge(input_tuple(), hash_accessor());
    m_output_tuple = return_data.m_output_tuple;

    return true;
}

std::vector<Tuple> TetEdgeCollapse::modified_primitives(PrimitiveType type) const
{
    if (type == PrimitiveType::Face) {
        return modified_triangles();
    } else {
        return {};
    }
}

std::string TetEdgeCollapse::name() const
{
    return "tet_mesh_collapse_edge";
}

Tuple TetEdgeCollapse::return_tuple() const
{
    return m_output_tuple;
}

std::vector<Tuple> TetEdgeCollapse::modified_triangles() const
{
    Simplex v(PrimitiveType::Vertex, m_output_tuple);
    auto sc = SimplicialComplex::open_star(mesh(), v);
    auto faces = sc.get_simplices(PrimitiveType::Face);
    std::vector<Tuple> ret;
    for (const auto& face : faces) {
        ret.emplace_back(face.tuple());
    }
    return ret;
}
} // namespace tet_mesh
} // namespace wmtk::operations
