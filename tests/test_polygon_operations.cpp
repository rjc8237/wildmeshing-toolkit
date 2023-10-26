#include <catch2/catch_test_macros.hpp>

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

    DEBUG_PolygonMesh m = digon();

    REQUIRE(m.is_connectivity_valid());
    OperationSettings<polygon_mesh::Splice> op_settings;
    bool success;

    SECTION("no change")
    {
        Tuple h = m.tuple_from_id(PH, 0);
        Tuple g = m.tuple_from_id(PH, 0);
        polygon_mesh::Splice op(m, h, g, op_settings);
        success = op();

        auto hn_accessor = m.create_base_accessor<long>(m.hn_handle());
        CHECK(hn_accessor.scalar_attribute(0) == 2);
        CHECK(hn_accessor.scalar_attribute(1) == 3);
        CHECK(hn_accessor.scalar_attribute(2) == 0);
        CHECK(hn_accessor.scalar_attribute(3) == 1);
        CHECK(m.get_all(PrimitiveType::Vertex).size() == 2);
        CHECK(m.get_all(PrimitiveType::Face).size() == 2);
    }
    SECTION("merge face merge vertex")
    {
        Tuple h = m.tuple_from_id(PH, 0);
        Tuple g = m.tuple_from_id(PH, 1);
        polygon_mesh::Splice op(m, h, g, op_settings);
        success = op();

        auto hn_accessor = m.create_base_accessor<long>(m.hn_handle());
        CHECK(hn_accessor.scalar_attribute(0) == 3);
        CHECK(hn_accessor.scalar_attribute(1) == 2);
        CHECK(hn_accessor.scalar_attribute(2) == 0);
        CHECK(hn_accessor.scalar_attribute(3) == 1);
        CHECK(m.get_all(PrimitiveType::Vertex).size() == 1);
        CHECK(m.get_all(PrimitiveType::Face).size() == 1);
    }
    SECTION("split face merge vertex")
    {
        Tuple h = m.tuple_from_id(PH, 0);
        Tuple g = m.tuple_from_id(PH, 2);
        polygon_mesh::Splice op(m, h, g, op_settings);
        success = op();

        auto hn_accessor = m.create_base_accessor<long>(m.hn_handle());
        CHECK(hn_accessor.scalar_attribute(0) == 0);
        CHECK(hn_accessor.scalar_attribute(1) == 3);
        CHECK(hn_accessor.scalar_attribute(2) == 2);
        CHECK(hn_accessor.scalar_attribute(3) == 1);
        CHECK(m.get_all(PrimitiveType::Vertex).size() == 1);
        CHECK(m.get_all(PrimitiveType::Face).size() == 3);
    }
    SECTION("merge face split vertex")
    {
        Tuple h = m.tuple_from_id(PH, 0);
        Tuple g = m.tuple_from_id(PH, 3);
        polygon_mesh::Splice op(m, h, g, op_settings);
        success = op();

        auto hn_accessor = m.create_base_accessor<long>(m.hn_handle());
        CHECK(hn_accessor.scalar_attribute(0) == 1);
        CHECK(hn_accessor.scalar_attribute(1) == 3);
        CHECK(hn_accessor.scalar_attribute(2) == 0);
        CHECK(hn_accessor.scalar_attribute(3) == 2);
        CHECK(m.get_all(PrimitiveType::Vertex).size() == 3);
        CHECK(m.get_all(PrimitiveType::Face).size() == 1);
    }

    CHECK(success);
    CHECK(m.get_all(PrimitiveType::Edge).size() == 2);
    CHECK(m.get_all(PrimitiveType::HalfEdge).size() == 4);
    CHECK(m.is_valid_connectivity());
}
