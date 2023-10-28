#include <stdlib.h>
#include <catch2/catch_test_macros.hpp>
#include <iostream>
#include <wmtk/PolygonMesh.hpp>
#include <wmtk/utils/Logger.hpp>
#include "tools/DEBUG_PolygonMesh.hpp"
#include "tools/PolygonMesh_examples.hpp"

using namespace wmtk;
using namespace wmtk::tests;

// TODO Make sure we have full coverage of methods

// Count the number of faces that are interior (i.e., not holes) in the mesh
long count_interior_faces(const DEBUG_PolygonMesh& m, const std::vector<Tuple> faces)
{
    long num_interior_faces = 0;
    for (const auto& face : faces) {
        if (!m.is_hole_face(face)) {
            ++num_interior_faces;
        }
    }
    return num_interior_faces;
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

// Compare vertex, edge, face, and halfedge indices of the tuples
bool tuple_equal(const DEBUG_PolygonMesh& m, const Tuple& t0, const Tuple& t1)
{
    const auto l = wmtk::logger().level();
    wmtk::logger().set_level(spdlog::level::err);
    const long v0 = m.id(t0, PrimitiveType::Vertex);
    const long e0 = m.id(t0, PrimitiveType::Edge);
    const long f0 = m.id(t0, PrimitiveType::Face);
    const long h0 = m.id(t0, PrimitiveType::HalfEdge);
    const long v1 = m.id(t1, PrimitiveType::Vertex);
    const long e1 = m.id(t1, PrimitiveType::Edge);
    const long f1 = m.id(t1, PrimitiveType::Face);
    const long h1 = m.id(t1, PrimitiveType::HalfEdge);
    wmtk::logger().set_level(l);
    return (v0 == v1) && (e0 == e1) && (f0 == f1) && (h0 == h1);
}

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
    // checking for every tuple t:
    // (1) t.switch_vertex().switch_vertex() == t
    // (2) t.switch_edge().switch_edge() == t
    // (3) t.switch_tri().switch_tri() == t

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

    size_t n_boundary_edges = 0;
    for (const Tuple& e : m.get_all(PrimitiveType::Edge)) {
        if (m.is_boundary(e)) {
            ++n_boundary_edges;
        }
    }
    CHECK(n_boundary_edges == 8);

    size_t n_boundary_vertices = 0;
    for (const Tuple& v : m.get_all(PrimitiveType::Vertex)) {
        if (m.is_boundary_vertex(v)) {
            ++n_boundary_vertices;
        }
    }
    CHECK(n_boundary_vertices == 8);

    const Tuple t1 = m.halfedge_tuple_from_vertex_in_face(0, 0);
    CHECK(m.is_boundary_edge(t1));
    CHECK(m.is_boundary_vertex(t1));

    const Tuple t2 = m.halfedge_tuple_from_vertex_in_face(1, 0);
    CHECK(!m.is_boundary_edge(t2));
    CHECK(m.is_boundary_vertex(t2));

    const Tuple t3 = m.halfedge_tuple_from_vertex_in_face(4, 0);
    CHECK(!m.is_boundary_edge(t3));
    CHECK(!m.is_boundary_vertex(t3));
}

// TODO Test get tuple
