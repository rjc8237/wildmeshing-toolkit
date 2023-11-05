#pragma once
#include <wmtk/PolygonMesh.hpp>
#include "PolygonMeshOperation.hpp"

namespace wmtk::operations::polygon_mesh {

/**
 * @class Atomic operations for changing connectivity
 */
class AtomicOperation : public PolygonMeshOperation
{
public:
    // constructor for default factory pattern construction
    AtomicOperation(Mesh& m);
    AtomicOperation(PolygonMesh& m);

protected:
    long new_face(bool is_hole);
    long new_edge();
    long new_vertex();

    void delete_face(long face_id);
    void delete_edge(long edge_id);
    void delete_vertex(long vertex_id);

    void set_face(long halfedge_id, long face_id);
    void set_next(long halfedge_id, long next_halfedge_id);
    void set_vertex(long halfedge_id, long vertex_id);

    long get_face(long halfedge_id);
    long get_next(long halfedge_id);
    long get_vertex(long halfedge_id);
    long get_halfedge_from_tuple(const Tuple& t);


private:
    Accessor<long> m_next_accessor;
    Accessor<long> m_prev_accessor;

    Accessor<long> m_to_accessor;
    Accessor<long> m_out_accessor;

    Accessor<long> m_hf_accessor;
    Accessor<long> m_fh_accessor;
    Accessor<char> m_f_is_hole_accessor;

    Accessor<char> m_v_flag_accessor;
    Accessor<char> m_e_flag_accessor;
    Accessor<char> m_f_flag_accessor;
    Accessor<char> m_h_flag_accessor;
};

} // namespace wmtk::operations::polygon_mesh
