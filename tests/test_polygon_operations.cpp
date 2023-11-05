#include <catch2/catch_test_macros.hpp>

#include <wmtk/operations/polygon_mesh/DeleteBubble.hpp>
#include <wmtk/operations/polygon_mesh/MakeEdge.hpp>
#include <wmtk/operations/polygon_mesh/Splice.hpp>
#include "tools/DEBUG_PolygonMesh.hpp"
#include "tools/PolygonMesh_examples.hpp"

using namespace wmtk;
using namespace wmtk::tests;

using PM = PolygonMesh;

constexpr PrimitiveType PV = PrimitiveType::Vertex;
constexpr PrimitiveType PE = PrimitiveType::Edge;
constexpr PrimitiveType PF = PrimitiveType::Face;
constexpr PrimitiveType PH = PrimitiveType::HalfEdge;

TEST_CASE("splice_edge_operation", "[operations][splice][polygon]")
{
    using namespace operations;

    DEBUG_PolygonMesh m;
    OperationSettings<polygon_mesh::Splice> op_settings;
    bool success;

    SECTION("no change")
    {
        m = digon();
        Tuple h = m.tuple_from_id(PH, 0);
        Tuple g = m.tuple_from_id(PH, 0);
        polygon_mesh::Splice op(m, h, g, op_settings);
        success = op();

        auto hn_accessor = m.create_base_accessor<long>(m.next_handle());
        CHECK(hn_accessor.scalar_attribute(0) == 2);
        CHECK(hn_accessor.scalar_attribute(1) == 3);
        CHECK(hn_accessor.scalar_attribute(2) == 0);
        CHECK(hn_accessor.scalar_attribute(3) == 1);
        CHECK(m.get_all(PrimitiveType::Vertex).size() == 2);
        CHECK(m.get_all(PrimitiveType::Face).size() == 2);
    }
    SECTION("merge face merge vertex")
    {
        m = digon();
        Tuple h = m.tuple_from_id(PH, 0);
        Tuple g = m.tuple_from_id(PH, 1);
        polygon_mesh::Splice op(m, h, g, op_settings);
        success = op();

        auto hn_accessor = m.create_base_accessor<long>(m.next_handle());
        CHECK(hn_accessor.scalar_attribute(0) == 3);
        CHECK(hn_accessor.scalar_attribute(1) == 2);
        CHECK(hn_accessor.scalar_attribute(2) == 0);
        CHECK(hn_accessor.scalar_attribute(3) == 1);
        CHECK(m.get_all(PrimitiveType::Vertex).size() == 1);
        CHECK(m.get_all(PrimitiveType::Face).size() == 1);
    }
    SECTION("split face merge vertex")
    {
        m = digon();
        Tuple h = m.tuple_from_id(PH, 0);
        Tuple g = m.tuple_from_id(PH, 2);
        polygon_mesh::Splice op(m, h, g, op_settings);
        success = op();

        auto hn_accessor = m.create_base_accessor<long>(m.next_handle());
        CHECK(hn_accessor.scalar_attribute(0) == 0);
        CHECK(hn_accessor.scalar_attribute(1) == 3);
        CHECK(hn_accessor.scalar_attribute(2) == 2);
        CHECK(hn_accessor.scalar_attribute(3) == 1);
        CHECK(m.get_all(PrimitiveType::Vertex).size() == 1);
        CHECK(m.get_all(PrimitiveType::Face).size() == 3);
    }
    SECTION("merge face split vertex")
    {
        m = digon();
        Tuple h = m.tuple_from_id(PH, 0);
        Tuple g = m.tuple_from_id(PH, 3);
        polygon_mesh::Splice op(m, h, g, op_settings);
        success = op();

        auto hn_accessor = m.create_base_accessor<long>(m.next_handle());
        CHECK(hn_accessor.scalar_attribute(0) == 1);
        CHECK(hn_accessor.scalar_attribute(1) == 3);
        CHECK(hn_accessor.scalar_attribute(2) == 0);
        CHECK(hn_accessor.scalar_attribute(3) == 2);
        CHECK(m.get_all(PrimitiveType::Vertex).size() == 3);
        CHECK(m.get_all(PrimitiveType::Face).size() == 1);
    }
    SECTION("split face split vertex")
    {
        m = torus();
        Tuple h = m.tuple_from_id(PH, 0);
        Tuple g = m.tuple_from_id(PH, 1);
        polygon_mesh::Splice op(m, h, g, op_settings);
        success = op();

        auto hn_accessor = m.create_base_accessor<long>(m.next_handle());
        CHECK(hn_accessor.scalar_attribute(0) == 3);
        CHECK(hn_accessor.scalar_attribute(1) == 2);
        CHECK(hn_accessor.scalar_attribute(2) == 1);
        CHECK(hn_accessor.scalar_attribute(3) == 0);
        CHECK(m.get_all(PrimitiveType::Vertex).size() == 2);
        CHECK(m.get_all(PrimitiveType::Face).size() == 2);
    }

    CHECK(success);
    CHECK(m.get_all(PrimitiveType::Edge).size() == 2);
    CHECK(m.get_all(PrimitiveType::HalfEdge).size() == 4);
    REQUIRE(m.is_connectivity_valid());
}

TEST_CASE("random_splice_edge_operation", "[operations][splice][polygon]")
{
    using namespace operations;

    long n_edges = 100;
    long n_halfedges = 2 * n_edges;
    DEBUG_PolygonMesh m = random_polygon_mesh(n_edges);
    OperationSettings<polygon_mesh::Splice> op_settings;

    // Randomly splice edges
    long num_splices = 100;
    for (long i = 0; i < num_splices; ++i) {
        long euler_characteristic_before = m.euler_characteristic();
        long num_components_before = m.count_simply_connected_components();
        Tuple h = m.tuple_from_id(PH, rand() % n_halfedges);
        Tuple g = m.tuple_from_id(PH, rand() % n_halfedges);
        polygon_mesh::Splice op(m, h, g, op_settings);
        bool is_topology_preserving = op.is_topology_preserving();
        bool success = op();

        CHECK(success);
        REQUIRE(m.is_connectivity_valid());
        CHECK(m.get_all(PrimitiveType::Edge).size() == n_edges);
        CHECK(m.get_all(PrimitiveType::HalfEdge).size() == (2 * n_edges));
        if (is_topology_preserving) {
            CHECK(m.euler_characteristic() == euler_characteristic_before);
            CHECK(m.count_simply_connected_components() == num_components_before);
        }
    }
}

TEST_CASE("make_edge_operation", "[operations][make_edge][polygon]")
{
    using namespace operations;

    DEBUG_PolygonMesh m;
    OperationSettings<polygon_mesh::MakeEdge> op_settings;
    bool success;

    m = digon();

    // Add first bubble
    polygon_mesh::MakeEdge op1(m, op_settings);
    success = op1();

    CHECK(success);
    REQUIRE(m.is_connectivity_valid());
    CHECK(m.euler_characteristic() == 4);
    CHECK(m.count_simply_connected_components() == 2);
    CHECK(m.get_all(PrimitiveType::Vertex).size() == 4);
    CHECK(m.get_all(PrimitiveType::Face).size() == 3);
    CHECK(m.get_all(PrimitiveType::Edge).size() == 3);
    CHECK(m.get_all(PrimitiveType::HalfEdge).size() == 6);

    // Add second bubble
    polygon_mesh::MakeEdge op2(m, op_settings);
    success = op2();

    CHECK(success);
    REQUIRE(m.is_connectivity_valid());
    CHECK(m.euler_characteristic() == 6);
    CHECK(m.count_simply_connected_components() == 3);
    CHECK(m.get_all(PrimitiveType::Vertex).size() == 6);
    CHECK(m.get_all(PrimitiveType::Face).size() == 4);
    CHECK(m.get_all(PrimitiveType::Edge).size() == 4);
    CHECK(m.get_all(PrimitiveType::HalfEdge).size() == 8);
}

TEST_CASE("delete_bubble_operation", "[operations][delete_bubble][polygon]")
{
    using namespace operations;

    DEBUG_PolygonMesh m;
    OperationSettings<polygon_mesh::DeleteBubble> op_settings;
    bool success;

    m = bubble();
    Tuple h = m.tuple_from_id(PrimitiveType::Edge, 0);

    // Delete bubble
    polygon_mesh::DeleteBubble op(m, h, op_settings);
    success = op();

    CHECK(success);
    REQUIRE(m.is_connectivity_valid());
    CHECK(m.euler_characteristic() == 0);
    CHECK(m.count_simply_connected_components() == 0);
    CHECK(m.get_all(PrimitiveType::Vertex).size() == 0);
    CHECK(m.get_all(PrimitiveType::Face).size() == 0);
    CHECK(m.get_all(PrimitiveType::Edge).size() == 0);
    CHECK(m.get_all(PrimitiveType::HalfEdge).size() == 0);
}

TEST_CASE("polygon_compound_atomic_operations", "[operations][polygon]")
{
    using namespace operations;

    DEBUG_PolygonMesh m;
    bool success;

    // Initialize one bubble
    m = bubble();
    Tuple h = m.tuple_from_id(PrimitiveType::Edge, 0);

    // Add another bubble
    OperationSettings<polygon_mesh::MakeEdge> op1_settings;
    polygon_mesh::MakeEdge op1(m, op1_settings);
    success = op1();
    Tuple g = op1.return_tuple();
    CHECK(success);
    REQUIRE(m.is_connectivity_valid());
    CHECK(m.euler_characteristic() == 4);
    CHECK(m.count_simply_connected_components() == 2);
    CHECK(m.get_all(PrimitiveType::Vertex).size() == 4);
    CHECK(m.get_all(PrimitiveType::Face).size() == 2);
    CHECK(m.get_all(PrimitiveType::Edge).size() == 2);
    CHECK(m.get_all(PrimitiveType::HalfEdge).size() == 4);

    // Splice first and second bubble
    OperationSettings<polygon_mesh::Splice> op2_settings;
    polygon_mesh::Splice op2(m, h, g, op2_settings);
    success = op2();
    CHECK(success);
    REQUIRE(m.is_connectivity_valid());
    CHECK(m.euler_characteristic() == 2);
    CHECK(m.count_simply_connected_components() == 1);
    CHECK(m.get_all(PrimitiveType::Vertex).size() == 3);
    CHECK(m.get_all(PrimitiveType::Face).size() == 1);
    CHECK(m.get_all(PrimitiveType::Edge).size() == 2);
    CHECK(m.get_all(PrimitiveType::HalfEdge).size() == 4);

    // TODO Might want to change tuples after splices
    // Undo splice
    OperationSettings<polygon_mesh::Splice> op3_settings;
    polygon_mesh::Splice op3(m, h, g, op3_settings);
    success = op3();
    CHECK(success);
    REQUIRE(m.is_connectivity_valid());
    CHECK(m.euler_characteristic() == 4);
    CHECK(m.count_simply_connected_components() == 2);
    CHECK(m.get_all(PrimitiveType::Vertex).size() == 4);
    CHECK(m.get_all(PrimitiveType::Face).size() == 2);
    CHECK(m.get_all(PrimitiveType::Edge).size() == 2);
    CHECK(m.get_all(PrimitiveType::HalfEdge).size() == 4);

    // Remove bubble
    OperationSettings<polygon_mesh::DeleteBubble> op4_settings;
    polygon_mesh::DeleteBubble op4(m, h, op4_settings);
    success = op4();
    CHECK(success);
    REQUIRE(m.is_connectivity_valid());
    CHECK(m.euler_characteristic() == 2);
    CHECK(m.count_simply_connected_components() == 1);
    CHECK(m.get_all(PrimitiveType::Vertex).size() == 2);
    CHECK(m.get_all(PrimitiveType::Face).size() == 1);
    CHECK(m.get_all(PrimitiveType::Edge).size() == 1);
    CHECK(m.get_all(PrimitiveType::HalfEdge).size() == 2);
}

// TODO Check hole face cases and preconditions