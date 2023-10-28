#pragma once

#include "Mesh.hpp"
#include "Tuple.hpp"

#include <Eigen/Core>

#include "Mesh.hpp"
#include "Tuple.hpp"

namespace wmtk {

class PolygonMesh : public Mesh
{
public:
    PolygonMesh();
    PolygonMesh(const PolygonMesh& o);
    PolygonMesh(PolygonMesh&& o);
    PolygonMesh& operator=(const PolygonMesh& o);
    PolygonMesh& operator=(PolygonMesh&& o);

    long top_cell_dimension() const override { return 2; }

    Tuple switch_tuple(const Tuple& tuple, PrimitiveType type) const override;

    /**
     * @brief Jump to the next halfedge by performing a switch of vertex and edge
     */
    Tuple next_halfedge(const Tuple& h_tuple) const;

    /**
     * @brief Jump to the previous halfedge by performing a switch of edge and vertex
     */
    Tuple prev_halfedge(const Tuple& h_tuple) const;

    /**
     * @brief Jump to the opposite halfedge by performing a switch of face and vertex
     */
    Tuple opp_halfedge(const Tuple& h_tuple) const;

    bool is_connectivity_valid() const override;
    bool is_local_connectivity_valid(const Tuple& v_tuple, PrimitiveType type) const;

    Tuple tuple_from_id(const PrimitiveType type, const long gid) const override;

    bool is_ccw(const Tuple& tuple) const override;

    /**
     * @brief Check if the face is a hole (i.e., not part of the surface)
     */
    bool is_hole_face(const Tuple& f_tuple) const;

    bool is_boundary(const Tuple& tuple) const override;
    bool is_boundary_vertex(const Tuple& tuple) const override;
    bool is_boundary_edge(const Tuple& tuple) const override;


    bool is_valid(const Tuple& tuple, ConstAccessor<long>& hash_accessor) const override;

    void initialize(
        Eigen::Ref<const VectorXl> next,
        Eigen::Ref<const VectorXl> prev,
        Eigen::Ref<const VectorXl> to,
        Eigen::Ref<const VectorXl> out,
        Eigen::Ref<const VectorXl> he2f,
        Eigen::Ref<const VectorXl> f2he);
    void initialize(Eigen::Ref<const VectorXl> next);
    void initialize_fv(const std::vector<std::vector<long>>& F);
    void initialize_fv(Eigen::Ref<const RowVectors3l> F);

protected:
    long id(const Tuple& tuple, PrimitiveType type) const override;


protected:
    attribute::MeshAttributeHandle<long> m_next_handle; // HalfEdge -> next HalfEdge
    attribute::MeshAttributeHandle<long> m_prev_handle; // HalfEdge -> previous HalfEdge

    attribute::MeshAttributeHandle<long> m_to_handle; // HalfEdge -> Vertex at tip
    attribute::MeshAttributeHandle<long> m_out_handle; // Vertex -> any outgoing HalfEdge

    attribute::MeshAttributeHandle<long> m_hf_handle; // HalfEdge -> adjacent Face
    attribute::MeshAttributeHandle<long> m_fh_handle; // Face -> any adjacent HalfEdge

    attribute::MeshAttributeHandle<char> m_f_is_hole_handle; // 1 if Face is a hole

    bool is_vertex_connectivity_valid(long vid) const;
    bool is_edge_connectivity_valid(long eid) const;
    bool is_face_connectivity_valid(long fid) const;
    bool is_halfedge_connectivity_valid(long hid) const;
};

} // namespace wmtk
