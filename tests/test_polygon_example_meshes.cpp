
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
    // the number of simplices (vertex, edge, then face)
    std::array<long, 3> simplex_counts = std::array<long, 3>{{-1, -1, -1}};
};

void run_debug_polygon_mesh(const DEBUG_PolygonMesh& m, const MeshDebugInfo& info)
{
    REQUIRE(m.is_connectivity_valid());
}
} // namespace

TEST_CASE("test_debug_polygon_meshes_triangle", "[examples]")
{
    DEBUG_PolygonMesh m = triangle();
    MeshDebugInfo info;
    run_debug_polygon_mesh(m, info);
}

TEST_CASE("test_debug_polygon_meshes_square", "[examples]")
{
    DEBUG_PolygonMesh m = square();
    MeshDebugInfo info;
    run_debug_polygon_mesh(m, info);
}

TEST_CASE("test_debug_polygon_meshes_pentagon", "[examples]")
{
    DEBUG_PolygonMesh m = pentagon();
    MeshDebugInfo info;
    run_debug_polygon_mesh(m, info);
}

TEST_CASE("test_debug_polygon_meshes_digon", "[examples]")
{
    DEBUG_PolygonMesh m = digon();
    MeshDebugInfo info;
    run_debug_polygon_mesh(m, info);
}

TEST_CASE("test_debug_polygon_meshes_monogon", "[examples]")
{
    DEBUG_PolygonMesh m = monogon();
    MeshDebugInfo info;
    run_debug_polygon_mesh(m, info);
}

TEST_CASE("test_debug_polygon_meshes_bubble", "[examples]")
{
    DEBUG_PolygonMesh m = bubble();
    MeshDebugInfo info;
    run_debug_polygon_mesh(m, info);
}

TEST_CASE("test_debug_polygon_meshes_dual_bubble", "[examples]")
{
    DEBUG_PolygonMesh m = dual_bubble();
    MeshDebugInfo info;
    run_debug_polygon_mesh(m, info);
}

TEST_CASE("test_debug_polygon_meshes_glued_polygons", "[examples]")
{
    DEBUG_PolygonMesh m = glued_polygons();
    MeshDebugInfo info;
    run_debug_polygon_mesh(m, info);
}

TEST_CASE("test_debug_polygon_meshes_grid", "[examples]")
{
    DEBUG_PolygonMesh m = grid();
    MeshDebugInfo info;
    run_debug_polygon_mesh(m, info);
}