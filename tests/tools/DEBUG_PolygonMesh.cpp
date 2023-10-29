#include "DEBUG_PolygonMesh.hpp"
#include <queue>
#include <stdexcept>
#include <wmtk/utils/Logger.hpp>

namespace wmtk::tests {

DEBUG_PolygonMesh::DEBUG_PolygonMesh(const PolygonMesh& m)
    : PolygonMesh(m)
{}
DEBUG_PolygonMesh::DEBUG_PolygonMesh(PolygonMesh&& m)
    : PolygonMesh(std::move(m))
{}


bool DEBUG_PolygonMesh::operator==(const DEBUG_PolygonMesh& o) const
{
    return static_cast<const PolygonMesh&>(*this) == static_cast<const PolygonMesh&>(o);
}
bool DEBUG_PolygonMesh::operator!=(const DEBUG_PolygonMesh& o) const
{
    return !(*this == o);
}

long DEBUG_PolygonMesh::id(const Tuple& tuple, PrimitiveType type) const
{
    return PolygonMesh::id(tuple, type);
}

Tuple DEBUG_PolygonMesh::halfedge_tuple_from_vertex_in_face(long vid, long fid) const
{
    Tuple f_tuple = tuple_from_id(PrimitiveType::Face, fid);
    Tuple f_tuple_iter = f_tuple;
    do {
        if (id(f_tuple_iter, PrimitiveType::Vertex) == vid) {
            break;
        }
        f_tuple_iter = next_halfedge(f_tuple_iter);
    } while (id(f_tuple_iter, PrimitiveType::HalfEdge) != id(f_tuple, PrimitiveType::HalfEdge));

    assert(id(f_tuple_iter, PrimitiveType::Vertex) == vid);
    return f_tuple_iter;
}

std::array<long, 4> DEBUG_PolygonMesh::primitive_counts() const
{
    return std::array<long, 4>{
        {long(get_all(PrimitiveType::Vertex).size()),
         long(get_all(PrimitiveType::Edge).size()),
         long(get_all(PrimitiveType::Face).size()),
         long(get_all(PrimitiveType::HalfEdge).size())}};
}

long DEBUG_PolygonMesh::find_simply_connected_components(std::vector<long>& component_ids) const
{
    component_ids.resize(capacity(PrimitiveType::HalfEdge), -1);
    const auto all_halfedge_tuples = get_all(PrimitiveType::HalfEdge);
    long component_id = 0;

    // BFS
    for (auto halfedge_tuple : all_halfedge_tuples) {
        const long hid = id(halfedge_tuple, PrimitiveType::HalfEdge);
        if (component_ids[hid] != -1) {
            continue;
        } // visited

        std::queue<long> q;
        q.push(hid);
        while (!q.empty()) {
            long cur_hid = q.front();
            q.pop();

            component_ids[cur_hid] = component_id;
            auto h_tuple = tuple_from_id(PrimitiveType::HalfEdge, cur_hid);

            // Add opposite halfedge to component
            auto opp_tuple = opp_halfedge(h_tuple);
            long opp_hid = id(opp_tuple, PrimitiveType::HalfEdge);
            if (component_ids[opp_hid] == -1) {
                q.push(opp_hid);
            }

            // Add next halfedge to component
            auto next_tuple = next_halfedge(h_tuple);
            long next_hid = id(next_tuple, PrimitiveType::HalfEdge);
            if (component_ids[next_hid] == -1) {
                q.push(next_hid);
            }
        }

        component_id++;
    }

    return component_id;
}

long DEBUG_PolygonMesh::count_simply_connected_components() const
{
    std::vector<long> component_ids;
    return find_simply_connected_components(component_ids);
}

long DEBUG_PolygonMesh::count_boundary_loops() const
{
    std::vector<long> boundary_loop_ids(capacity(PrimitiveType::Edge), -1);
    const auto all_edge_tuples = get_all(PrimitiveType::Edge);
    long boundary_loop_id = 0;
    for (auto edge_tuple : all_edge_tuples) {
        if (!is_boundary(edge_tuple)) {
            continue;
        }
        const long eid = id(edge_tuple, PrimitiveType::Edge);
        if (boundary_loop_ids[eid] != -1) {
            continue;
        } // visited

        Tuple cur_edge = edge_tuple;
        long cur_eid = id(cur_edge, PrimitiveType::Edge);
        // find one boundary loop
        while (boundary_loop_ids[eid] == -1 || cur_eid != eid) {
            boundary_loop_ids[cur_eid] = boundary_loop_id;

            // find next boundary edge
            cur_edge = switch_edge(switch_vertex(cur_edge));
            while (!is_boundary(cur_edge)) {
                cur_edge = switch_edge(switch_face(cur_edge));
            }

            cur_eid = id(cur_edge, PrimitiveType::Edge);
        }

        boundary_loop_id++;
    }

    return boundary_loop_id;
}

long DEBUG_PolygonMesh::count_hole_faces() const
{
    auto faces = get_all(PrimitiveType::Face);
    long num_hole_faces = 0;
    for (const auto& face : faces) {
        if (is_hole_face(face)) {
            ++num_hole_faces;
        }
    }

    return num_hole_faces;
}

long DEBUG_PolygonMesh::euler_characteristic() const
{
    auto counts = primitive_counts();
    long num_boundary_loops = count_boundary_loops();
    long num_hole_faces = count_hole_faces();
    return counts[0] - counts[1] + counts[2] + num_boundary_loops - num_hole_faces;
}

long DEBUG_PolygonMesh::genus() const
{
    if (count_simply_connected_components() != 1) {
        wmtk::logger().error("not a connected surface\n");
        throw std::runtime_error("GenusComputeError");
    }
    return (2 - euler_characteristic()) / 2;
}

} // namespace wmtk::tests
