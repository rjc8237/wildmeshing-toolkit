#pragma once
#include <array>
#include <wmtk/PolygonMesh.hpp>

namespace wmtk::tests {
class DEBUG_PolygonMesh : public PolygonMesh //, public virtual DEBUG_Mesh
{
public:
    using PolygonMesh::PolygonMesh;
    DEBUG_PolygonMesh(const PolygonMesh& m);
    DEBUG_PolygonMesh(PolygonMesh&& m);
    using PolygonMesh::operator=;


    bool operator==(const DEBUG_PolygonMesh& o) const;
    bool operator!=(const DEBUG_PolygonMesh& o) const;

    Tuple halfedge_tuple_from_vertex_in_face(long vid, long fid) const;


    long id(const Tuple& tuple, PrimitiveType type) const override;
    using PolygonMesh::tuple_from_id;

    std::array<long, 4> primitive_counts() const;
    long find_simply_connected_components(std::vector<long>& component_ids) const;
    long count_simply_connected_components() const;
    long count_boundary_loops() const;
    long count_hole_faces() const;
    long euler_characteristic() const;
    long genus() const;
};

} // namespace wmtk::tests
