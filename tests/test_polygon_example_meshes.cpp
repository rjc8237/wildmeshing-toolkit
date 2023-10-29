
#include <catch2/catch_test_macros.hpp>
#include "tools/DEBUG_PolygonMesh.hpp"
#include "tools/PolygonMesh_examples.hpp"

using namespace wmtk;
using namespace wmtk::tests;

// TODO Add additional checks like elements counts and genus
// It would be better to make the computation accessible for other tests like splice

namespace {
struct MeshDebugInfo
{
    std::string name = "RENAME ME";
    long boundary_curves = -1; // number of distinct boundary curves in the mesh
    long genus = -1; // TODO (also what definition of genus?)
    long simply_connected_components = -1; // TODO (face-face connected topologies)
    bool is_oriented = false; // TODO (make sure nface neighbors use opposite orientation
    // the number of simplices (vertex, edge, face, then halfedge)
    std::array<long, 4> primitive_counts = std::array<long, 4>{{-1, -1, -1, -1}};
};

void run_debug_polygon_mesh(const DEBUG_PolygonMesh& m, const MeshDebugInfo& info)
{
    // Validate info
    REQUIRE(info.boundary_curves >= 0);
    REQUIRE(info.simply_connected_components >= 0);
    for (const long count : info.primitive_counts) {
        REQUIRE(count >= 0);
    }

    // Check topology
    REQUIRE(m.is_connectivity_valid());
    if (info.genus >= 0) {
        CHECK(m.genus() == info.genus);
    }
    CHECK(m.count_simply_connected_components() == info.simply_connected_components);
    CHECK(m.primitive_counts() == info.primitive_counts);
    CHECK(m.count_boundary_loops() == info.boundary_curves);
}
} // namespace

TEST_CASE("test_debug_polygon_meshes_triangle", "[examples]")
{
    DEBUG_PolygonMesh m = triangle();
    MeshDebugInfo info;
    info.name = "triangle";
    info.genus = 0;
    info.boundary_curves = 0;
    info.simply_connected_components = 1;
    info.primitive_counts = std::array<long, 4>{{3, 3, 2, 6}};
    run_debug_polygon_mesh(m, info);
}

TEST_CASE("test_debug_polygon_meshes_square", "[examples]")
{
    DEBUG_PolygonMesh m = square();
    MeshDebugInfo info;
    info.name = "square";
    info.genus = 0;
    info.boundary_curves = 0;
    info.simply_connected_components = 1;
    info.primitive_counts = std::array<long, 4>{{4, 4, 2, 8}};
    run_debug_polygon_mesh(m, info);
}

TEST_CASE("test_debug_polygon_meshes_pentagon", "[examples]")
{
    DEBUG_PolygonMesh m = pentagon();
    MeshDebugInfo info;
    info.name = "pentagon";
    info.genus = 0;
    info.boundary_curves = 0;
    info.simply_connected_components = 1;
    info.primitive_counts = std::array<long, 4>{{5, 5, 2, 10}};
    run_debug_polygon_mesh(m, info);
}

TEST_CASE("test_debug_polygon_meshes_digon", "[examples]")
{
    DEBUG_PolygonMesh m = digon();
    MeshDebugInfo info;
    info.name = "digon";
    info.genus = 0;
    info.boundary_curves = 0;
    info.simply_connected_components = 1;
    info.primitive_counts = std::array<long, 4>{{2, 2, 2, 4}};
    run_debug_polygon_mesh(m, info);
}

TEST_CASE("test_debug_polygon_meshes_monogon", "[examples]")
{
    DEBUG_PolygonMesh m = monogon();
    MeshDebugInfo info;
    info.name = "monogon";
    info.genus = 0;
    info.boundary_curves = 0;
    info.simply_connected_components = 1;
    info.primitive_counts = std::array<long, 4>{{1, 1, 2, 2}};
    run_debug_polygon_mesh(m, info);
}

TEST_CASE("test_debug_polygon_meshes_bubble", "[examples]")
{
    DEBUG_PolygonMesh m = bubble();
    MeshDebugInfo info;
    info.name = "bubble";
    info.genus = 0;
    info.boundary_curves = 0;
    info.simply_connected_components = 1;
    info.primitive_counts = std::array<long, 4>{{2, 1, 1, 2}};
    run_debug_polygon_mesh(m, info);
}

TEST_CASE("test_debug_polygon_meshes_dual_bubble", "[examples]")
{
    DEBUG_PolygonMesh m = dual_bubble();
    MeshDebugInfo info;
    info.name = "dual_bubble";
    info.genus = 0;
    info.boundary_curves = 0;
    info.simply_connected_components = 1;
    info.primitive_counts = std::array<long, 4>{{1, 1, 2, 2}};
    run_debug_polygon_mesh(m, info);
}

TEST_CASE("test_debug_polygon_meshes_glued_polygons", "[examples]")
{
    DEBUG_PolygonMesh m = glued_polygons();
    MeshDebugInfo info;
    info.name = "glued_polygons";
    info.genus = 0;
    info.boundary_curves = 1;
    info.simply_connected_components = 1;
    info.primitive_counts = std::array<long, 4>{{6, 7, 3, 14}};
    run_debug_polygon_mesh(m, info);
}

TEST_CASE("test_debug_polygon_meshes_torus", "[examples]")
{
    DEBUG_PolygonMesh m = torus();
    MeshDebugInfo info;
    info.name = "torus";
    info.genus = 1;
    info.boundary_curves = 0;
    info.simply_connected_components = 1;
    info.primitive_counts = std::array<long, 4>{{1, 2, 1, 4}};
    run_debug_polygon_mesh(m, info);
}

TEST_CASE("test_debug_polygon_meshes_annulus", "[examples]")
{
    DEBUG_PolygonMesh m = annulus();
    MeshDebugInfo info;
    info.name = "annulus";
    info.genus = 0;
    info.boundary_curves = 2;
    info.simply_connected_components = 1;
    info.primitive_counts = std::array<long, 4>{{8, 9, 3, 18}};
    run_debug_polygon_mesh(m, info);
}

TEST_CASE("test_debug_polygon_meshes_two_squares", "[examples]")
{
    DEBUG_PolygonMesh m = two_squares();
    MeshDebugInfo info;
    info.name = "two_squares";
    info.boundary_curves = 2;
    info.simply_connected_components = 2;
    info.primitive_counts = std::array<long, 4>{{8, 8, 4, 16}};
    run_debug_polygon_mesh(m, info);
}

TEST_CASE("test_debug_polygon_meshes_grid", "[examples]")
{
    DEBUG_PolygonMesh m = grid();
    MeshDebugInfo info;
    info.name = "triangle";
    info.genus = 0;
    info.boundary_curves = 1;
    info.simply_connected_components = 1;
    info.primitive_counts = std::array<long, 4>{{9, 12, 5, 24}};
    run_debug_polygon_mesh(m, info);
}