#pragma once

#include "Mesh.hpp"
#include "Tuple.hpp"

#include <Eigen/Core>
#include "Mesh.hpp"

#include "Mesh.hpp"
#include "Tuple.hpp"

namespace wmtk {

namespace operations::polygon_mesh {
class AtomicOperation;
} // namespace operations::polygon_mesh

/**
 * @brief Representation of a general oriented manifold polygonal mesh.
 *
 * The geometric primitives of the polygonal mesh are vertices, edges, faces, and halfedges.
 *
 */
class PolygonMesh : public Mesh
{
public:
    friend class operations::polygon_mesh::AtomicOperation;

    PolygonMesh();
    PolygonMesh(const PolygonMesh& o);
    PolygonMesh(PolygonMesh&& o);
    PolygonMesh& operator=(const PolygonMesh& o);
    PolygonMesh& operator=(PolygonMesh&& o);

    long top_cell_dimension() const override { return 2; }

    Tuple switch_tuple(const Tuple& tuple, PrimitiveType type) const override;
    Tuple tuple_from_id(const PrimitiveType type, const long gid) const override;

    bool is_connectivity_valid() const override;
    bool is_valid(const Tuple& tuple, ConstAccessor<long>& hash_accessor) const override;

    bool is_ccw(const Tuple& tuple) const override;

    bool is_boundary(const Tuple& tuple) const override;
    bool is_boundary_vertex(const Tuple& tuple) const override;
    bool is_boundary_edge(const Tuple& tuple) const override;

    /**
     * @brief Jump to the next halfedge by performing a switch of vertex and edge
     *
     * @param h_tuple: tuple specifying a halfedge
     * @return tuple corresponding to the next halfedge
     */
    Tuple next_halfedge(const Tuple& h_tuple) const;

    /**
     * @brief Jump to the previous halfedge by performing a switch of edge and vertex
     *
     * @param h_tuple: tuple specifying a halfedge
     * @return tuple corresponding to the previous halfedge
     */
    Tuple prev_halfedge(const Tuple& h_tuple) const;

    /**
     * @brief Jump to the opposite halfedge by performing a switch of face and vertex
     *
     * @param h_tuple: tuple specifying a halfedge
     * @return tuple corresponding to the opposite halfedge
     */
    Tuple opp_halfedge(const Tuple& h_tuple) const;

    /**
     * @brief Check if a face is a hole (i.e., not part of the surface)
     *
     * @param f_tuple: tuple specifying a face
     * @return true iff the face is a hole
     */
    bool is_hole_face(const Tuple& f_tuple) const;

    /**
     * @brief Initialize a polygon mesh from a full halfedge representation with implicit opp.
     *
     * @param next: size #he vector, next halfedge id
     * @param prev: size #he vector, prev halfedge id
     * @param to: size #he vector, halfedge vertex tip id
     * @param out: size #v vector, arbitrary halfedge id outgoing from vertex
     * @param he2f: size #he vector, face id adjacent to halfedge
     * @param f2he: size #f vector, arbitrary halfedge id adjacent to face
     */
    void initialize(
        Eigen::Ref<const VectorXl> next,
        Eigen::Ref<const VectorXl> prev,
        Eigen::Ref<const VectorXl> to,
        Eigen::Ref<const VectorXl> out,
        Eigen::Ref<const VectorXl> he2f,
        Eigen::Ref<const VectorXl> f2he);

    /**
     * @brief Initialize a polygon mesh from a minimal halfedge representation with implicit opp.
     *
     * Face and vertex ids are generated implicitly as needed.
     *
     * @param next: size #he vector, next halfedge id
     */
    void initialize(Eigen::Ref<const VectorXl> next);

    /**
     * @brief Initialize a polygon mesh from an FV representation.
     *
     * @param F: list of lists of face vertices
     */
    void initialize_fv(const std::vector<std::vector<long>>& F);

    /**
     * @brief Initialize a triangle mesh from an FV representation.
     *
     * @param F: $fx3 matrix of face vertices
     */
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
};

} // namespace wmtk
