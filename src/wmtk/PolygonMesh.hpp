#pragma once

#include "Mesh.hpp"
#include "Tuple.hpp"

#include <Eigen/Core>

#include "Mesh.hpp"
#include "Tuple.hpp"

/**
 * Proposed Implementation of the General Polygonal Mesh Data Structure
 *
 * A general polygonal mesh is represented using tuples of a consistently oriented face, edge,
 * and vertex together with orientation-reversing switch vertex, switch edge, and switch face maps.
 *
 * Using the tuple structure, the global cid is a halfedge (equivalently, an edge with a normal
 * direction), and a local vid. The vertex index can be interpreted as the tip of the edge and is
 * thus equivalent to specifying a vertex of the edge as well as a direction, which together with
 * the normal direction specifies an orientation of the face. The local eid and local fid are both
 * set to default -1.
 *
 * The implementation is done using the standard halfedge operations, next and previous halfedge,
 * encoded as pairs of edge attributes (two per edge, one for each halfedge). Opp is inferred from
 * the pairing of halfedges into edges. The switch vertex, switch edge, and switch face operations
 * can then be implemented as simple combinations of next, prev, and opp along with the additional
 * operation of switching the vertex of the halfedge. Switching the vertex corresponds to switching
 * the orientation of the surface, so next becomes prev and vice versa.
 *
 * Note that halfedge only supports oriented surfaces. However, the oriented halfedge could easily
 * be extended to support nonorientable surfaces by adding an additional edge attributes to indicate
 * whether the two adjacent faces should be glued with consistent or anti-consistent orientation.
 *
 * The choice of halfedges directed with the local vertex index is natural for the given tuple
 * interface as the halfedge is a 1-simplex in the given face, and using halfedge allows for some
 * theoretical simplifications as well as more compact memory usage.
 *
 * Also note that the natural generalization of this implementation is the facet-edge data
 * structure of Dobkin and Laszlo, where the primitive is an oriented face-edge
 * triple. However, while the halfedge can be encoded as an edge attribute, the facet-edge
 * data needs to be encoded on actual pairs of faces and edges as each face can be incident
 * to many edges and vice versa.
 */

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
};

} // namespace wmtk
