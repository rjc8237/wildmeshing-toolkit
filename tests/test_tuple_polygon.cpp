#include <stdlib.h>
#include <catch2/catch_test_macros.hpp>
#include <iostream>
#include <wmtk/PolygonMesh.hpp>
#include <wmtk/utils/Logger.hpp>
#include "tools/DEBUG_PolygonMesh.hpp"
#include "tools/PolygonMesh_examples.hpp"

using namespace wmtk;
using namespace wmtk::tests;

// Count the number of faces that are interior (i.e., not holes) in the mesh
long count_interior_faces(const DEBUG_PolygonMesh& m, const std::vector<Tuple> faces)
{
    return m.get_all(PrimitiveType::Face).size() - m.count_hole_faces();
}


TEST_CASE("polygon_initialize", "[mesh_creation],[tuple_polygon]")
{
    DEBUG_PolygonMesh m = triangle();

    // Check connectivity is valid
    REQUIRE(m.is_connectivity_valid());

    // Get all primitives for the mesh
    const std::vector<Tuple> vertices = m.get_all(PrimitiveType::Vertex);
    CHECK(vertices.size() == 3);
    const std::vector<Tuple> edges = m.get_all(PrimitiveType::Edge);
    CHECK(edges.size() == 3);
    const std::vector<Tuple> faces = m.get_all(PrimitiveType::Face);
    CHECK(count_interior_faces(m, faces) == 2);
    const std::vector<Tuple> halfedges = m.get_all(PrimitiveType::HalfEdge);
    CHECK(halfedges.size() == 6);

    // Check tuple validity
    CHECK(m.is_valid_slow(vertices[0]));
    CHECK(m.is_valid_slow(edges[0]));
    CHECK(m.is_valid_slow(faces[0]));
    CHECK(m.is_valid_slow(halfedges[0]));
}

// Check the number of primitive types, and check that tuple_from_id and id compose to the identity
void check_tuple_generation(const DEBUG_PolygonMesh& m, long n_vertices, long n_edges, long n_faces)
{
    SECTION("vertices")
    {
        const std::vector<Tuple> vertices = m.get_all(PrimitiveType::Vertex);
        REQUIRE(vertices.size() == n_vertices);
        for (long vi = 0; vi < n_vertices; ++vi) {
            CHECK(m.id(vertices[vi], PrimitiveType::Vertex) == vi);
        }
    }
    SECTION("edges")
    {
        const std::vector<Tuple> edges = m.get_all(PrimitiveType::Edge);
        REQUIRE(edges.size() == n_edges);
        for (long ei = 0; ei < n_edges; ++ei) {
            CHECK(m.id(edges[ei], PrimitiveType::Edge) == ei);
        }
    }
    SECTION("faces")
    {
        const std::vector<Tuple> faces = m.get_all(PrimitiveType::Face);
        CHECK(count_interior_faces(m, faces) == n_faces);
        for (long fi = 0; fi < n_faces; ++fi) {
            CHECK(m.id(faces[fi], PrimitiveType::Face) == fi);
        }
    }
    SECTION("halfedges")
    {
        const std::vector<Tuple> halfedges = m.get_all(PrimitiveType::HalfEdge);
        REQUIRE(halfedges.size() == 2 * n_edges);
        for (long hi = 0; hi < 2 * n_edges; ++hi) {
            CHECK(m.id(halfedges[hi], PrimitiveType::HalfEdge) == hi);
        }
    }
}

TEST_CASE("polygon_1_triangle", "[tuple_generation],[tuple_polygon]")
{
    DEBUG_PolygonMesh m = triangle();
    check_tuple_generation(m, 3, 3, 2);
}

TEST_CASE("glued_polygons", "[tuple_generation],[tuple_polygon]")
{
    DEBUG_PolygonMesh m = glued_polygons();
    check_tuple_generation(m, 6, 7, 2);
}

// Randomly switch a vertex, edge, or face of a tuple
void random_switch(const DEBUG_PolygonMesh& m, Tuple& t)
{
    switch (rand() % 3) {
    case 0: t = m.switch_tuple(t, PrimitiveType::Vertex); break;
    case 1: t = m.switch_tuple(t, PrimitiveType::Edge); break;
    case 2: t = m.switch_tuple(t, PrimitiveType::Face); break;
    default: break;
    }
}


TEST_CASE("polygon_random_switches", "[tuple_operation],[tuple_polygon]")
{
    // Check that random tuple switches are valid on a random mesh

    DEBUG_PolygonMesh m = random_polygon_mesh(20);

    SECTION("vertices")
    {
        const std::vector<Tuple> vertex_tuples = m.get_all(PrimitiveType::Vertex);
        for (size_t i = 0; i < vertex_tuples.size(); ++i) {
            Tuple t = vertex_tuples[i];
            for (size_t j = 0; j < 10; j++) {
                random_switch(m, t);
                CHECK(m.is_valid_slow(t));
            }
        }
    }

    SECTION("edges")
    {
        const std::vector<Tuple> edge_tuples = m.get_all(PrimitiveType::Edge);
        for (size_t i = 0; i < edge_tuples.size(); ++i) {
            Tuple t = edge_tuples[i];
            for (size_t j = 0; j < 10; j++) {
                random_switch(m, t);
                CHECK(m.is_valid_slow(t));
            }
        }
    }

    SECTION("faces")
    {
        const std::vector<Tuple> face_tuples = m.get_all(PrimitiveType::Face);
        for (size_t i = 0; i < face_tuples.size(); ++i) {
            Tuple t = face_tuples[i];
            for (size_t j = 0; j < 10; j++) {
                random_switch(m, t);
                CHECK(m.is_valid_slow(t));
            }
        }
    }
}

TEST_CASE("polygon_known_switches", "[tuple_operation],[tuple_polygon]")
{
    // Check explicit traversal of a mesh with known face and vertex ids

    DEBUG_PolygonMesh m = grid();
    Tuple t = m.halfedge_tuple_from_vertex_in_face(0, 0);
    CHECK(m.id(t, PrimitiveType::Vertex) == 0);
    CHECK(m.id(t, PrimitiveType::Face) == 0);
    bool initial_is_ccw = m.is_ccw(t);

    t = m.switch_tuple(t, PrimitiveType::Vertex);
    CHECK(m.id(t, PrimitiveType::Vertex) == 1);
    CHECK(m.id(t, PrimitiveType::Face) == 0);
    CHECK(m.is_ccw(t) != initial_is_ccw);

    t = m.switch_tuple(t, PrimitiveType::Edge);
    CHECK(m.id(t, PrimitiveType::Vertex) == 1);
    CHECK(m.id(t, PrimitiveType::Face) == 0);
    CHECK(m.is_ccw(t) == initial_is_ccw);

    t = m.switch_tuple(t, PrimitiveType::Face);
    CHECK(m.id(t, PrimitiveType::Vertex) == 1);
    CHECK(m.id(t, PrimitiveType::Face) == 1);
    CHECK(m.is_ccw(t) != initial_is_ccw);

    t = m.next_halfedge(t);
    CHECK(m.id(t, PrimitiveType::Vertex) == 4);
    CHECK(m.id(t, PrimitiveType::Face) == 1);
    CHECK(m.is_ccw(t) != initial_is_ccw);

    t = m.switch_tuple(t, PrimitiveType::Vertex);
    CHECK(m.id(t, PrimitiveType::Vertex) == 5);
    CHECK(m.id(t, PrimitiveType::Face) == 1);
    CHECK(m.is_ccw(t) == initial_is_ccw);

    t = m.prev_halfedge(t);
    CHECK(m.id(t, PrimitiveType::Vertex) == 2);
    CHECK(m.id(t, PrimitiveType::Face) == 1);
    CHECK(m.is_ccw(t) == initial_is_ccw);

    t = m.opp_halfedge(t);
    CHECK(m.id(t, PrimitiveType::Vertex) == 5);
    CHECK(m.id(t, PrimitiveType::Face) == 4);
    CHECK(m.is_ccw(t) == initial_is_ccw);
}


// Compare vertex, edge, face, and halfedge indices of the tuples
bool tuple_equal(const DEBUG_PolygonMesh& m, const Tuple& t0, const Tuple& t1)
{
    const long v0 = m.id(t0, PrimitiveType::Vertex);
    const long e0 = m.id(t0, PrimitiveType::Edge);
    const long f0 = m.id(t0, PrimitiveType::Face);
    const long h0 = m.id(t0, PrimitiveType::HalfEdge);

    const long v1 = m.id(t1, PrimitiveType::Vertex);
    const long e1 = m.id(t1, PrimitiveType::Edge);
    const long f1 = m.id(t1, PrimitiveType::Face);
    const long h1 = m.id(t1, PrimitiveType::HalfEdge);

    return (v0 == v1) && (e0 == e1) && (f0 == f1) && (h0 == h1);
}

TEST_CASE("polygon_fixed_point_switches", "[tuple_operation],[tuple_polygon]")
{
    // Check that tuple switches are valid and do not have fixed points on meshes with a single
    // vertex, edge, and face

    DEBUG_PolygonMesh m;

    SECTION("vertex_fixed")
    {
        m = dual_bubble();
        const std::vector<Tuple> vertex_tuples = m.get_all(PrimitiveType::Vertex);
        for (const auto& t : vertex_tuples) {
            Tuple t_switch = m.switch_tuple(t, PrimitiveType::Vertex);
            CHECK(m.is_valid_slow(t_switch));
            CHECK(!t.same_ids(t_switch));
            CHECK(m.id(t_switch, PrimitiveType::Vertex) == m.id(t, PrimitiveType::Vertex));
        }
    }

    SECTION("edge_fixed")
    {
        m = bubble();
        const std::vector<Tuple> edge_tuples = m.get_all(PrimitiveType::Edge);
        for (const auto& t : edge_tuples) {
            Tuple t_switch = m.switch_tuple(t, PrimitiveType::Edge);
            CHECK(m.is_valid_slow(t_switch));
            CHECK(!t.same_ids(t_switch));
            CHECK(m.id(t_switch, PrimitiveType::Edge) == m.id(t, PrimitiveType::Edge));
        }
    }

    SECTION("face_fixed")
    {
        m = bubble();
        const std::vector<Tuple> face_tuples = m.get_all(PrimitiveType::Face);
        for (const auto& t : face_tuples) {
            Tuple t_switch = m.switch_tuple(t, PrimitiveType::Face);
            CHECK(m.is_valid_slow(t_switch));
            CHECK(!t.same_ids(t_switch));
            CHECK(m.id(t_switch, PrimitiveType::Face) == m.id(t, PrimitiveType::Face));
        }
    }
}

// Check that primitive switches are their own inverse for a given tuple
void check_double_switches(const DEBUG_PolygonMesh& m, const Tuple& t)
{
    const Tuple t_after_v = m.switch_vertex(m.switch_vertex(t));
    CHECK(tuple_equal(m, t, t_after_v));
    const Tuple t_after_e = m.switch_edge(m.switch_edge(t));
    CHECK(tuple_equal(m, t, t_after_e));
    const Tuple t_after_f = m.switch_face(m.switch_face(t));
    CHECK(tuple_equal(m, t, t_after_f));
}

TEST_CASE("polygon_double_switches", "[tuple_operation],[tuple_polygon]")
{
    DEBUG_PolygonMesh m = random_polygon_mesh(20);

    SECTION("vertices")
    {
        const std::vector<Tuple> vertices = m.get_all(PrimitiveType::Vertex);
        for (const auto& t : vertices) {
            check_double_switches(m, t);
        }
    }
    SECTION("edges")
    {
        const std::vector<Tuple> edges = m.get_all(PrimitiveType::Edge);
        for (const auto& t : edges) {
            check_double_switches(m, t);
        }
    }
    SECTION("faces")
    {
        const std::vector<Tuple> faces = m.get_all(PrimitiveType::Face);
        for (const auto& t : faces) {
            check_double_switches(m, t);
        }
    }
}

// Check that next and prev are inverse for a given tuple
void check_next_prev(const DEBUG_PolygonMesh& m, const Tuple& t)
{
    const Tuple t_pn = m.next_halfedge(m.prev_halfedge(t));
    CHECK(tuple_equal(m, t, t_pn));
    const Tuple t_np = m.prev_halfedge(m.next_halfedge(t));
    CHECK(tuple_equal(m, t, t_np));
}

TEST_CASE("polygon_next_prev", "[tuple_operation],[tuple_polygon]")
{
    DEBUG_PolygonMesh m = random_polygon_mesh(50);

    SECTION("vertices")
    {
        const std::vector<Tuple> vertices = m.get_all(PrimitiveType::Vertex);
        for (const Tuple& t : vertices) {
            check_next_prev(m, t);
        }
    }
    SECTION("edges")
    {
        const std::vector<Tuple> edges = m.get_all(PrimitiveType::Edge);
        for (const Tuple& t : edges) {
            check_next_prev(m, t);
        }
    }
    SECTION("faces")
    {
        const std::vector<Tuple> faces = m.get_all(PrimitiveType::Face);
        for (const Tuple& t : faces) {
            check_next_prev(m, t);
        }
    }
}

// Check that opp is its own inverse for a given tuple
void check_opp_opp(const DEBUG_PolygonMesh& m, const Tuple& t)
{
    const Tuple t_iter = m.opp_halfedge(m.opp_halfedge(t));
    CHECK(tuple_equal(m, t, t_iter));
}

TEST_CASE("polygon_opp_opp", "[tuple_operation],[tuple_polygon]")
{
    DEBUG_PolygonMesh m = random_polygon_mesh(50);

    SECTION("vertices")
    {
        const std::vector<Tuple> vertices = m.get_all(PrimitiveType::Vertex);
        for (const Tuple& t : vertices) {
            check_opp_opp(m, t);
        }
    }
    SECTION("edges")
    {
        const std::vector<Tuple> edges = m.get_all(PrimitiveType::Edge);
        for (const Tuple& t : edges) {
            check_opp_opp(m, t);
        }
    }
    SECTION("faces")
    {
        const std::vector<Tuple> faces = m.get_all(PrimitiveType::Face);
        for (const Tuple& t : faces) {
            check_opp_opp(m, t);
        }
    }
}

// Check that next-next-next is the identity for a given tuple (for a triangle)
void check_next_next_next(const DEBUG_PolygonMesh& m, const Tuple& t)
{
    const Tuple t_iter = m.next_halfedge(m.next_halfedge(m.next_halfedge(t)));
    CHECK(tuple_equal(m, t, t_iter));
}

TEST_CASE("polygon_next_next_next", "[tuple_operation],[tuple_polygon]")
{
    DEBUG_PolygonMesh m = triangle();

    SECTION("vertices")
    {
        const std::vector<Tuple> vertices = m.get_all(PrimitiveType::Vertex);
        for (const Tuple& t : vertices) {
            check_next_next_next(m, t);
        }
    }
    SECTION("edges")
    {
        const std::vector<Tuple> edges = m.get_all(PrimitiveType::Edge);
        for (const Tuple& t : edges) {
            check_next_next_next(m, t);
        }
    }
    SECTION("faces")
    {
        const std::vector<Tuple> faces = m.get_all(PrimitiveType::Face);
        for (const Tuple& t : faces) {
            check_next_next_next(m, t);
        }
    }
}

TEST_CASE("polygon_one_ring_iteration", "[tuple_operation],[tuple_polygon]")
{
    DEBUG_PolygonMesh m;
    {
        RowVectors3l tris;
        tris.resize(6, 3);
        tris.row(0) = Eigen::Matrix<long, 3, 1>{0, 1, 2};
        tris.row(1) = Eigen::Matrix<long, 3, 1>{0, 2, 3};
        tris.row(2) = Eigen::Matrix<long, 3, 1>{0, 3, 4};
        tris.row(3) = Eigen::Matrix<long, 3, 1>{0, 4, 5};
        tris.row(4) = Eigen::Matrix<long, 3, 1>{0, 5, 6};
        tris.row(5) = Eigen::Matrix<long, 3, 1>{0, 6, 1};
        m.initialize_fv(tris);
    }

    // Check that alternately switching edges and faces circulates around a vertex
    const std::vector<Tuple> vertices = m.get_all(PrimitiveType::Vertex);
    REQUIRE(vertices.size() == 7);
    for (const auto& t : vertices) {
        // face-edge switch
        Tuple t_iter = t;
        for (size_t i = 0; i < 6; ++i) {
            t_iter = m.switch_tuple(t_iter, PrimitiveType::Face);
            t_iter = m.switch_tuple(t_iter, PrimitiveType::Edge);
            if (tuple_equal(m, t, t_iter)) {
                break;
            }
        }
        CHECK(tuple_equal(m, t, t_iter));

        // edge-face switch
        t_iter = t;
        for (size_t i = 0; i < 6; ++i) {
            t_iter = m.switch_tuple(t_iter, PrimitiveType::Edge);
            t_iter = m.switch_tuple(t_iter, PrimitiveType::Face);
            if (tuple_equal(m, t, t_iter)) {
                break;
            }
        }
        CHECK(tuple_equal(m, t, t_iter));
    }
}

TEST_CASE("polygon_is_boundary", "[tuple_polygon]")
{
    DEBUG_PolygonMesh m = grid();

    // Check boundary edge count
    size_t n_boundary_edges = 0;
    for (const Tuple& e : m.get_all(PrimitiveType::Edge)) {
        if (m.is_boundary(e)) {
            ++n_boundary_edges;
        }
    }
    CHECK(n_boundary_edges == 8);

    // Check boundary vertex count
    size_t n_boundary_vertices = 0;
    for (const Tuple& v : m.get_all(PrimitiveType::Vertex)) {
        if (m.is_boundary_vertex(v)) {
            ++n_boundary_vertices;
        }
    }
    CHECK(n_boundary_vertices == 8);

    // Check a tuple with vertex and edge on the boundary
    const Tuple t1 = m.halfedge_tuple_from_vertex_in_face(0, 0);
    CHECK(m.is_boundary_edge(t1));
    CHECK(m.is_boundary_vertex(t1));

    // Check a tuple with vertex but not edge on the boundary
    const Tuple t2 = m.halfedge_tuple_from_vertex_in_face(1, 0);
    CHECK(!m.is_boundary_edge(t2));
    CHECK(m.is_boundary_vertex(t2));

    // Check a tuple with neither vertex nor edge on the boundary
    const Tuple t3 = m.halfedge_tuple_from_vertex_in_face(4, 0);
    CHECK(!m.is_boundary_edge(t3));
    CHECK(!m.is_boundary_vertex(t3));

    // Check hole face
    const Tuple t4 = m.opp_halfedge(t1);
    CHECK(!m.is_hole_face(t1));
    CHECK(m.is_hole_face(t4));
}

// TODO Test get tuple
