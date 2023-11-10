#include <catch2/catch_test_macros.hpp>

#include <wmtk/operations/polygon_mesh/DeleteBubble.hpp>
#include <wmtk/operations/polygon_mesh/FillHole.hpp>
#include <wmtk/operations/polygon_mesh/JoinFace.hpp>
#include <wmtk/operations/polygon_mesh/JoinVertex.hpp>
#include <wmtk/operations/polygon_mesh/MakeEdge.hpp>
#include <wmtk/operations/polygon_mesh/MakeHole.hpp>
#include <wmtk/operations/polygon_mesh/Splice.hpp>
#include <wmtk/operations/polygon_mesh/SplitFace.hpp>
#include <wmtk/operations/polygon_mesh/SplitVertex.hpp>
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

TEST_CASE("make_hole_operations", "[operations][make_hole][polygon]")
{
    using namespace operations;
    DEBUG_PolygonMesh m = grid();
    bool success;

    // Make a valid hole face
    Tuple f0 = m.tuple_from_id(PF, 0);
    OperationSettings<polygon_mesh::MakeHole> op1_settings;
    polygon_mesh::MakeHole op1(m, f0, op1_settings);
    CHECK(op1.precondition());
    success = op1();

    CHECK(success);
    REQUIRE(m.is_connectivity_valid());
    CHECK(m.count_simply_connected_components() == 1);
    CHECK(m.genus() == 0);
    CHECK(m.count_boundary_loops() == 1);
    CHECK(m.count_hole_faces() == 2);
    CHECK(m.count_interior_faces() == 3);

    // Make a invalid hole face
    Tuple f3 = m.tuple_from_id(PF, 3);
    OperationSettings<polygon_mesh::MakeHole> op2_settings;
    polygon_mesh::MakeHole op2(m, f3, op2_settings);
    CHECK(!op2.precondition());
}

TEST_CASE("fill_hole_operations", "[operations][fill_hole][polygon]")
{
    using namespace operations;
    DEBUG_PolygonMesh m = grid();
    bool success;

    // Make all faces holes
    for (long i = 0; i < 4; ++i) {
        Tuple fi = m.tuple_from_id(PF, i);
        OperationSettings<polygon_mesh::MakeHole> op_settings;
        polygon_mesh::MakeHole op(m, fi, op_settings);
        CHECK(op.precondition());
        success = op();
        CHECK(success);
        REQUIRE(m.is_connectivity_valid());
        CHECK(m.count_hole_faces() == 2 + i);
        CHECK(m.count_interior_faces() == 3 - i);
    }

    // Fill a valid hole
    Tuple f0 = m.tuple_from_id(PF, 0);
    OperationSettings<polygon_mesh::FillHole> make_hole_op_settings;
    polygon_mesh::FillHole op1(m, f0, make_hole_op_settings);
    CHECK(op1.precondition());
    success = op1();

    CHECK(success);
    REQUIRE(m.is_connectivity_valid());
    CHECK(m.count_simply_connected_components() == 1);
    CHECK(m.genus() == 0);
    CHECK(m.count_boundary_loops() == 1);
    CHECK(m.count_hole_faces() == 4);
    CHECK(m.count_interior_faces() == 1);

    // Fill a invalid hole face
    Tuple f3 = m.tuple_from_id(PF, 3);
    OperationSettings<polygon_mesh::FillHole> op2_settings;
    polygon_mesh::FillHole op2(m, f3, op2_settings);
    CHECK(!op2.precondition());
}


TEST_CASE("open_splice_edge_operation", "[operations][splice][polygon]")
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
        polygon_mesh::MakeHole(m, g, OperationSettings<polygon_mesh::MakeHole>())();
        polygon_mesh::Splice op(m, h, g, op_settings);
        assert(op.precondition());
        success = op();

        CHECK(success);
        CHECK(m.get_all(PrimitiveType::Vertex).size() == 2);
        CHECK(m.get_all(PrimitiveType::Face).size() == 2);
        CHECK(m.count_hole_faces() == 1);
    }
    SECTION("merge single first hole face")
    {
        m = digon();
        Tuple h = m.tuple_from_id(PH, 0);
        Tuple g = m.tuple_from_id(PH, 1);
        polygon_mesh::MakeHole(m, h, OperationSettings<polygon_mesh::MakeHole>())();
        polygon_mesh::Splice op(m, h, g, op_settings);
        assert(op.precondition());
        success = op();

        CHECK(success);
        CHECK(m.get_all(PrimitiveType::Vertex).size() == 1);
        CHECK(m.get_all(PrimitiveType::Face).size() == 1);
        CHECK(m.count_hole_faces() == 0);
    }
    SECTION("merge single second hole face")
    {
        m = digon();
        Tuple h = m.tuple_from_id(PH, 0);
        Tuple g = m.tuple_from_id(PH, 1);
        polygon_mesh::MakeHole(m, g, OperationSettings<polygon_mesh::MakeHole>())();
        polygon_mesh::Splice op(m, h, g, op_settings);
        assert(op.precondition());
        success = op();

        CHECK(success);
        CHECK(m.get_all(PrimitiveType::Vertex).size() == 1);
        CHECK(m.get_all(PrimitiveType::Face).size() == 1);
        CHECK(m.count_hole_faces() == 0);
    }
    SECTION("merge double hole face")
    {
        m = digon();
        Tuple h = m.tuple_from_id(PH, 0);
        Tuple g = m.tuple_from_id(PH, 1);
        polygon_mesh::MakeHole(m, h, OperationSettings<polygon_mesh::MakeHole>())();
        polygon_mesh::MakeHole(m, g, OperationSettings<polygon_mesh::MakeHole>())();
        polygon_mesh::Splice op(m, h, g, op_settings);
        assert(op.precondition());
        success = op();

        CHECK(success);
        CHECK(m.get_all(PrimitiveType::Vertex).size() == 1);
        CHECK(m.get_all(PrimitiveType::Face).size() == 1);
        CHECK(m.count_hole_faces() == 1);
    }
    SECTION("split hole face")
    {
        m = digon();
        Tuple h = m.tuple_from_id(PH, 0);
        Tuple g = m.tuple_from_id(PH, 2);
        Tuple f0 = m.tuple_from_id(PF, 0);
        Tuple f1 = m.tuple_from_id(PF, 1);
        polygon_mesh::MakeHole(m, f0, OperationSettings<polygon_mesh::MakeHole>())();
        polygon_mesh::MakeHole(m, f1, OperationSettings<polygon_mesh::MakeHole>())();
        polygon_mesh::Splice op(m, h, g, op_settings);
        assert(op.precondition());
        success = op();

        CHECK(success);
        CHECK(m.get_all(PrimitiveType::Vertex).size() == 1);
        CHECK(m.get_all(PrimitiveType::Face).size() == 3);
        CHECK(m.count_hole_faces() == 3);
    }
    SECTION("invalid split")
    {
        m = digon();
        Tuple h = m.tuple_from_id(PH, 0);
        Tuple g = m.tuple_from_id(PH, 2);
        polygon_mesh::MakeHole(m, h, OperationSettings<polygon_mesh::MakeHole>())();
        polygon_mesh::Splice op(m, h, g, op_settings);
        assert(!op.precondition());
    }

    CHECK(m.get_all(PrimitiveType::Edge).size() == 2);
    CHECK(m.get_all(PrimitiveType::HalfEdge).size() == 4);
    REQUIRE(m.is_connectivity_valid());
}

// Helper function to determine if sub is a contiguous subsequence of vec starting
// at some given offset into vec (with periodic indexing)
bool is_periodic_subsequence(
    const std::vector<long>& sub,
    const std::vector<long>& vec,
    long offset)
{
    long n = vec.size();
    long m = sub.size();
    if (n < m) {
        return false; // Cannot have longer subsequence than vector
    }
    if (m == 0) {
        return true; // Empty subsequence
    }

    for (long j = 0; j < m; ++j) {
        if (sub[j] != vec[(j + offset) % n]) {
            return false;
        }
    }

    return true;
}

// Helper function to determine if sub is a contiguous subsequence of vec starting
// at any offset into vec (with periodic indexing)
bool is_periodic_subsequence(const std::vector<long>& sub, const std::vector<long>& vec)
{
    // Try all possible starting points in the vector for the subsequence
    long n = vec.size();
    for (long i = 0; i < n; ++i) {
        if (is_periodic_subsequence(sub, vec, i)) {
            return true;
        }
    }

    return false;
}

// Check the postconditions for the split face operation
bool split_face_postcondition(
    const DEBUG_PolygonMesh& m0,
    const DEBUG_PolygonMesh& m,
    const Tuple& e_tuple,
    const Tuple& h_tuple,
    const Tuple& g_tuple)
{
    // Get input halfedges, their next halfedges in the input mesh, and new halfedge ids
    std::array<long, 2> hids = {m0.id(h_tuple, PH), m0.id(g_tuple, PH)};
    std::array<long, 2> hnids = {
        m0.id(m0.next_halfedge(h_tuple), PH),
        m0.id(m0.next_halfedge(g_tuple), PH)};
    std::array<long, 2> eids = {m.id(e_tuple, PH), m.id(m.opp_halfedge(e_tuple), PH)};

    // Check face loops
    for (long i = 0; i < 2; ++i) {
        std::vector<long> f0 = m0.get_face_id_loop_range(hnids[1 - i], hnids[i]);
        std::vector<long> f = m.get_face_id_loop(eids[i]);
        if ((f.size() != f0.size() + 1) || (!is_periodic_subsequence(f0, f, 1))) {
            return false;
        }
    }

    // Check distinct vertex loops
    if (m0.id(h_tuple, PV) != m0.id(g_tuple, PV)) {
        for (long i = 0; i < 2; ++i) {
            std::vector<long> v0 = m0.get_vertex_id_loop(hids[i]);
            std::vector<long> v = m.get_vertex_id_loop(eids[1 - i]);
            if ((v.size() != v0.size() + 1) || (!is_periodic_subsequence(v0, v, 1))) {
                return false;
            }
        }
    }
    // Check common vertex loops
    else {
        std::vector<long> v0 = m0.get_vertex_id_loop(hids[0]);
        std::vector<long> v = m.get_vertex_id_loop(eids[0]);
        if (v.size() != v0.size() + 2) {
            return false;
        }
        for (long i = 0; i < 2; ++i) {
            std::vector<long> v0 = m0.get_vertex_id_loop_range(hids[i], hids[1 - i]);
            std::vector<long> v = m.get_vertex_id_loop_range(eids[1 - i], eids[i]);
            if ((v.size() != v0.size() + 1) || (!is_periodic_subsequence(v0, v, 1))) {
                return false;
            }
        }
    }

    // Check the split faces are holes iff the original face is a hole
    bool is_hole = m0.is_hole_face(h_tuple);
    if ((m.is_hole_face(e_tuple) != is_hole) ||
        (m.is_hole_face(m.opp_halfedge(e_tuple)) != is_hole)) {
        return false;
    }

    return true;
}

TEST_CASE("split_face_operation", "[operations][split_face][polygon]")
{
    using namespace operations;

    DEBUG_PolygonMesh m0, m;
    OperationSettings<polygon_mesh::SplitFace> op_settings;
    bool success;

    SECTION("split digon")
    {
        m0 = digon();
        m = digon();
        Tuple h = m.tuple_from_id(PH, 0);
        Tuple g = m.tuple_from_id(PH, 2);
        polygon_mesh::SplitFace op(m, h, g, op_settings);
        CHECK(op.precondition());
        success = op();
        Tuple e = op.return_tuple();

        CHECK(split_face_postcondition(m0, m, e, h, g));
        CHECK(m.get_all(PrimitiveType::Vertex).size() == 2);
        CHECK(m.get_all(PrimitiveType::Edge).size() == 3);
        CHECK(m.get_all(PrimitiveType::Face).size() == 3);
        CHECK(m.count_hole_faces() == 0);
        CHECK(m.count_interior_faces() == 3);
    }
    SECTION("split triangle")
    {
        m0 = triangle();
        m = triangle();
        Tuple h = m.tuple_from_id(PH, 0);
        Tuple g = m.tuple_from_id(PH, 2);
        polygon_mesh::SplitFace op(m, h, g, op_settings);
        CHECK(op.precondition());
        success = op();
        Tuple e = op.return_tuple();

        CHECK(split_face_postcondition(m0, m, e, h, g));
        CHECK(m.get_all(PrimitiveType::Vertex).size() == 3);
        CHECK(m.get_all(PrimitiveType::Edge).size() == 4);
        CHECK(m.get_all(PrimitiveType::Face).size() == 3);
        CHECK(m.count_hole_faces() == 0);
        CHECK(m.count_interior_faces() == 3);
    }
    SECTION("split face with common vertex")
    {
        m0 = torus();
        m = torus();
        Tuple h = m.tuple_from_id(PH, 0);
        Tuple g = m.tuple_from_id(PH, 1);
        polygon_mesh::SplitFace op(m, h, g, op_settings);
        CHECK(op.precondition());
        success = op();
        Tuple e = op.return_tuple();

        CHECK(split_face_postcondition(m0, m, e, h, g)); // TODO Check this by hand
        CHECK(m.get_all(PrimitiveType::Vertex).size() == 1);
        CHECK(m.get_all(PrimitiveType::Edge).size() == 3);
        CHECK(m.get_all(PrimitiveType::Face).size() == 2);
        CHECK(m.count_hole_faces() == 0);
        CHECK(m.count_interior_faces() == 2);
    }
    SECTION("split hole face")
    {
        m0 = triangle();
        polygon_mesh::MakeHole(m0, m0.tuple_from_id(PH, 0), {})();
        m = triangle();
        polygon_mesh::MakeHole(m, m.tuple_from_id(PH, 0), {})();

        Tuple h = m.tuple_from_id(PH, 0);
        Tuple g = m.tuple_from_id(PH, 4);
        polygon_mesh::SplitFace op(m, h, g, op_settings);
        CHECK(op.precondition());
        success = op();
        Tuple e = op.return_tuple();

        CHECK(split_face_postcondition(m0, m, e, h, g));
        CHECK(m.get_all(PrimitiveType::Vertex).size() == 3);
        CHECK(m.get_all(PrimitiveType::Edge).size() == 4);
        CHECK(m.get_all(PrimitiveType::Face).size() == 3);
        CHECK(m.count_hole_faces() == 2);
        CHECK(m.count_interior_faces() == 1);
    }
    SECTION("invalid different face split")
    {
        m = triangle();

        Tuple h = m.tuple_from_id(PH, 0);
        Tuple g = m.tuple_from_id(PH, 1);
        polygon_mesh::SplitFace op(m, h, g, op_settings);
        CHECK(!op.precondition());
    }
    SECTION("invalid same halfedge split")
    {
        m = triangle();

        Tuple h = m.tuple_from_id(PH, 0);
        Tuple g = m.tuple_from_id(PH, 0);
        polygon_mesh::SplitFace op(m, h, g, op_settings);
        CHECK(!op.precondition());
    }

    CHECK(success);
    REQUIRE(m.is_connectivity_valid());
}

// Check the postconditions for the join face operation
bool join_face_postcondition(
    const DEBUG_PolygonMesh& m0,
    const DEBUG_PolygonMesh& m,
    const Tuple& h_tuple)
{
    std::array<long, 2> hids = {m0.id(h_tuple, PH), m0.id(m0.opp_halfedge(h_tuple), PH)};
    std::array<long, 2> hnids = {
        m0.id(m0.next_halfedge(h_tuple), PH),
        m0.id(m0.next_halfedge(m0.opp_halfedge(h_tuple)), PH)};
    std::array<long, 2> hpids = {
        m0.id(m0.prev_halfedge(h_tuple), PH),
        m0.id(m0.prev_halfedge(m0.opp_halfedge(h_tuple)), PH)};

    // monogon face removed by join
    // TODO Implement postcondition for this case if we allow this in the precondition
    if (hids[0] == hnids[0]) {
        return true;
    }
    if (hids[1] == hnids[1]) {
        return true;
    }

    // Check face loops
    for (long i = 0; i < 2; ++i) {
        std::vector<long> f0 = m0.get_face_id_loop(hids[i]);
        std::vector<long> f = m.get_face_id_loop_range(hnids[i], hnids[1 - i]);
        if ((f0.size() != f.size() + 1) || (!is_periodic_subsequence(f, f0, 1))) {
            return false;
        }
    }

    // Check distinct vertex loops
    if (m0.id(h_tuple, PV) != m0.id(m0.opp_halfedge(h_tuple), PV)) {
        for (long i = 0; i < 2; ++i) {
            std::vector<long> v0 = m0.get_vertex_id_loop(hids[1 - i]);
            std::vector<long> v = m.get_vertex_id_loop(hpids[i]);
            if ((v0.size() != v.size() + 1) || (!is_periodic_subsequence(v, v0, 1))) {
                return false;
            }
        }
    }
    // Check common vertex loops
    else {
        std::vector<long> v0 = m0.get_vertex_id_loop(hids[0]);
        std::vector<long> v = m.get_vertex_id_loop(hpids[0]);
        if (v0.size() != v.size() + 2) {
            return false;
        }
        for (long i = 0; i < 2; ++i) {
            std::vector<long> v0 = m0.get_vertex_id_loop_range(hids[i], hids[1 - i]);
            std::vector<long> v = m.get_vertex_id_loop_range(hpids[1 - i], hpids[i]);
            if ((v0.size() != v.size() + 1) || (!is_periodic_subsequence(v, v0, 1))) {
                return false;
            }
        }
    }

    // Check the joined face is a hole iff the original faces are both a hole
    bool is_hole = (m0.is_hole_face(h_tuple) && m0.is_hole_face(m0.opp_halfedge(h_tuple)));
    if ((m.is_hole_face(m.tuple_from_id(PH, hnids[0])) != is_hole) ||
        (m.is_hole_face(m.tuple_from_id(PH, hnids[1])) != is_hole)) {
        return false;
    }

    return true;
}

TEST_CASE("join_face_operation", "[operations][split_vertex][polygon]")
{
    using namespace operations;

    DEBUG_PolygonMesh m0, m;
    OperationSettings<polygon_mesh::JoinFace> op_settings;
    bool success;

    SECTION("join digon face")
    {
        m0 = digon();
        m = digon();

        Tuple h = m.tuple_from_id(PH, 0);
        polygon_mesh::JoinFace op(m, h, op_settings);
        CHECK(op.precondition());
        success = op();

        CHECK(join_face_postcondition(m0, m, h));
        CHECK(m.get_all(PrimitiveType::Vertex).size() == 2);
        CHECK(m.get_all(PrimitiveType::Edge).size() == 1);
        CHECK(m.get_all(PrimitiveType::Face).size() == 1);
        CHECK(m.count_hole_faces() == 0);
        CHECK(m.count_interior_faces() == 1);
    }
    SECTION("join triangle face")
    {
        m0 = triangle();
        m = triangle();

        Tuple h = m.tuple_from_id(PH, 2);
        polygon_mesh::JoinFace op(m, h, op_settings);
        CHECK(op.precondition());
        success = op();

        CHECK(join_face_postcondition(m0, m, h));
        CHECK(m.get_all(PrimitiveType::Vertex).size() == 3);
        CHECK(m.get_all(PrimitiveType::Edge).size() == 2);
        CHECK(m.get_all(PrimitiveType::Face).size() == 1);
        CHECK(m.count_hole_faces() == 0);
        CHECK(m.count_interior_faces() == 1);
    }
    SECTION("join leaf face")
    {
        m0 = monogon_with_leaf();
        m = monogon_with_leaf();

        Tuple h = m.tuple_from_id(PH, 0);
        polygon_mesh::JoinFace op(m, h, op_settings);
        CHECK(op.precondition());
        success = op();

        CHECK(join_face_postcondition(m0, m, h));
        CHECK(m.get_all(PrimitiveType::Vertex).size() == 2);
        CHECK(m.get_all(PrimitiveType::Edge).size() == 1);
        CHECK(m.get_all(PrimitiveType::Face).size() == 1);
        CHECK(m.count_hole_faces() == 0);
        CHECK(m.count_interior_faces() == 1);
    }
    SECTION("join face with hole face")
    {
        m0 = triangle();
        polygon_mesh::MakeHole(m0, m0.tuple_from_id(PH, 0), {})();
        m = triangle();
        polygon_mesh::MakeHole(m, m.tuple_from_id(PH, 0), {})();
        Tuple h = m.tuple_from_id(PH, 2);

        polygon_mesh::JoinFace op(m, h, op_settings);
        CHECK(op.precondition());
        success = op();

        CHECK(join_face_postcondition(m0, m, h));
        CHECK(m.get_all(PrimitiveType::Vertex).size() == 3);
        CHECK(m.get_all(PrimitiveType::Edge).size() == 2);
        CHECK(m.get_all(PrimitiveType::Face).size() == 1);
        CHECK(m.count_hole_faces() == 0);
        CHECK(m.count_interior_faces() == 1);
    }
    SECTION("join face with two hole face")
    {
        m0 = triangle();
        polygon_mesh::MakeHole(m0, m0.tuple_from_id(PH, 0), {})();
        polygon_mesh::MakeHole(m0, m0.tuple_from_id(PH, 1), {})();
        m = triangle();
        polygon_mesh::MakeHole(m, m.tuple_from_id(PH, 0), {})();
        polygon_mesh::MakeHole(m, m.tuple_from_id(PH, 1), {})();

        Tuple h = m.tuple_from_id(PH, 2);
        polygon_mesh::JoinFace op(m, h, op_settings);
        CHECK(op.precondition());
        success = op();

        CHECK(join_face_postcondition(m0, m, h));
        CHECK(m.get_all(PrimitiveType::Vertex).size() == 3);
        CHECK(m.get_all(PrimitiveType::Edge).size() == 2);
        CHECK(m.get_all(PrimitiveType::Face).size() == 1);
        CHECK(m.count_hole_faces() == 1);
        CHECK(m.count_interior_faces() == 0);
    }
    SECTION("invalid same face join")
    {
        m = torus();

        Tuple h = m.tuple_from_id(PH, 0);
        polygon_mesh::JoinFace op(m, h, op_settings);
        CHECK(!op.precondition());
    }

    CHECK(success);
    REQUIRE(m.is_connectivity_valid());
}

bool split_vertex_postcondition(
    const DEBUG_PolygonMesh& m0,
    const DEBUG_PolygonMesh& m,
    const Tuple& e_tuple,
    const Tuple& h_tuple,
    const Tuple& g_tuple)
{
    // Get input halfedge, their next halfedges in the input mesh, and new halfedge ids
    std::array<long, 2> hids = {m0.id(h_tuple, PH), m0.id(g_tuple, PH)};
    std::array<long, 2> hnids = {
        m0.id(m0.next_halfedge(h_tuple), PH),
        m0.id(m0.next_halfedge(g_tuple), PH)};
    std::array<long, 2> eids = {m.id(e_tuple, PH), m.id(m.opp_halfedge(e_tuple), PH)};

    // Check vertex loops
    for (long i = 0; i < 2; ++i) {
        std::vector<long> v0 = m0.get_vertex_id_loop_range(hids[i], hids[1 - i]);
        std::vector<long> v = m.get_vertex_id_loop(eids[1 - i]);
        if ((v.size() != v0.size() + 1) || (!is_periodic_subsequence(v0, v, 1))) {
            return false;
        }
    }

    // Check distinct face loops
    if (m0.id(h_tuple, PF) != m0.id(g_tuple, PF)) {
        for (long i = 0; i < 2; ++i) {
            std::vector<long> f0 = m0.get_face_id_loop(hnids[i]);
            std::vector<long> f = m.get_face_id_loop(eids[i]);
            if ((f.size() != f0.size() + 1) || (!is_periodic_subsequence(f0, f, 1))) {
                return false;
            }
        }
    }
    // Check common face loops
    else {
        std::vector<long> f0 = m0.get_face_id_loop(hnids[0]);
        std::vector<long> f = m.get_face_id_loop(eids[0]);
        if (f.size() != f0.size() + 2) {
            return false;
        }
        for (long i = 0; i < 2; ++i) {
            std::vector<long> f0 = m0.get_face_id_loop_range(hnids[i], hnids[1 - i]);
            std::vector<long> f = m.get_face_id_loop_range(eids[i], eids[1 - i]);
            if ((f.size() != f0.size() + 1) || (!is_periodic_subsequence(f0, f, 1))) {
                return false;
            }
        }
    }

    // Check the new faces are holes iff the original faces are holes
    if ((m.is_hole_face(e_tuple) != m0.is_hole_face(h_tuple)) ||
        (m.is_hole_face(m.opp_halfedge(e_tuple)) != m0.is_hole_face(g_tuple))) {
        return false;
    }

    return true;
}

TEST_CASE("split_vertex_operation", "[operations][split_vertex][polygon]")
{
    using namespace operations;

    DEBUG_PolygonMesh m0, m;
    OperationSettings<polygon_mesh::SplitVertex> op_settings;
    bool success;

    SECTION("split digon vertex")
    {
        m0 = digon();
        m = digon();
        Tuple h = m.tuple_from_id(PH, 0);
        Tuple g = m.tuple_from_id(PH, 3);
        polygon_mesh::SplitVertex op(m, h, g, op_settings);
        CHECK(op.precondition());
        success = op();
        Tuple e = op.return_tuple();

        CHECK(split_vertex_postcondition(m0, m, e, h, g));
        CHECK(m.get_all(PrimitiveType::Vertex).size() == 3);
        CHECK(m.get_all(PrimitiveType::Edge).size() == 3);
        CHECK(m.get_all(PrimitiveType::Face).size() == 2);
        CHECK(m.count_hole_faces() == 0);
        CHECK(m.count_interior_faces() == 2);
    }
    SECTION("split triangle vertex")
    {
        m0 = triangle();
        m = triangle();
        Tuple h = m.tuple_from_id(PH, 2);
        Tuple g = m.tuple_from_id(PH, 5);
        polygon_mesh::SplitVertex op(m, h, g, op_settings);
        CHECK(op.precondition());
        success = op();
        Tuple e = op.return_tuple();

        CHECK(split_vertex_postcondition(m0, m, e, h, g));
        CHECK(m.get_all(PrimitiveType::Vertex).size() == 4);
        CHECK(m.get_all(PrimitiveType::Edge).size() == 4);
        CHECK(m.get_all(PrimitiveType::Face).size() == 2);
        CHECK(m.count_hole_faces() == 0);
        CHECK(m.count_interior_faces() == 2);
    }
    SECTION("split vertex with common face")
    {
        m0 = torus();
        m = torus();
        Tuple h = m.tuple_from_id(PH, 0);
        Tuple g = m.tuple_from_id(PH, 1);
        // TODO Try g = 1,2,3
        // Tuple g = m.tuple_from_id(PH, 2);
        polygon_mesh::SplitVertex op(m, h, g, op_settings);
        CHECK(op.precondition());
        success = op();
        Tuple e = op.return_tuple();

        CHECK(split_vertex_postcondition(m0, m, e, h, g));
        CHECK(m.get_all(PrimitiveType::Vertex).size() == 2);
        CHECK(m.get_all(PrimitiveType::Edge).size() == 3);
        CHECK(m.get_all(PrimitiveType::Face).size() == 1);
        CHECK(m.count_hole_faces() == 0);
        CHECK(m.count_interior_faces() == 1);
    }
    SECTION("split vertex with hole face")
    {
        m0 = triangle();
        polygon_mesh::MakeHole(m0, m0.tuple_from_id(PH, 0), {})();
        m = triangle();
        polygon_mesh::MakeHole(m, m.tuple_from_id(PH, 0), {})();

        Tuple h = m.tuple_from_id(PH, 0);
        Tuple g = m.tuple_from_id(PH, 3);
        polygon_mesh::SplitVertex op(m, h, g, op_settings);
        CHECK(op.precondition());
        success = op();
        Tuple e = op.return_tuple();

        CHECK(split_vertex_postcondition(m0, m, e, h, g));
        CHECK(m.get_all(PrimitiveType::Vertex).size() == 4);
        CHECK(m.get_all(PrimitiveType::Edge).size() == 4);
        CHECK(m.get_all(PrimitiveType::Face).size() == 2);
        CHECK(m.count_hole_faces() == 1);
        CHECK(m.count_interior_faces() == 1);
    }
    SECTION("split vertex with two hole faces")
    {
        m0 = triangle();
        polygon_mesh::MakeHole(m0, m0.tuple_from_id(PH, 0), {})();
        polygon_mesh::MakeHole(m0, m0.tuple_from_id(PH, 1), {})();
        m = triangle();
        polygon_mesh::MakeHole(m, m.tuple_from_id(PH, 0), {})();
        polygon_mesh::MakeHole(m, m.tuple_from_id(PH, 1), {})();

        Tuple h = m.tuple_from_id(PH, 0);
        Tuple g = m.tuple_from_id(PH, 3);
        polygon_mesh::SplitVertex op(m, h, g, op_settings);
        CHECK(op.precondition());
        success = op();
        Tuple e = op.return_tuple();

        CHECK(split_vertex_postcondition(m0, m, e, h, g));
        CHECK(m.get_all(PrimitiveType::Vertex).size() == 4);
        CHECK(m.get_all(PrimitiveType::Edge).size() == 4);
        CHECK(m.get_all(PrimitiveType::Face).size() == 2);
        CHECK(m.count_hole_faces() == 2);
        CHECK(m.count_interior_faces() == 0);
    }
    SECTION("invalid different vertex split")
    {
        m = triangle();

        Tuple h = m.tuple_from_id(PH, 0);
        Tuple g = m.tuple_from_id(PH, 1);
        polygon_mesh::SplitVertex op(m, h, g, op_settings);
        CHECK(!op.precondition());
    }
    SECTION("invalid same halfedge split")
    {
        m = triangle();

        Tuple h = m.tuple_from_id(PH, 0);
        Tuple g = m.tuple_from_id(PH, 0);
        polygon_mesh::SplitVertex op(m, h, g, op_settings);
        CHECK(!op.precondition());
    }

    CHECK(success);
    REQUIRE(m.is_connectivity_valid());
}

bool join_vertex_postcondition(
    const DEBUG_PolygonMesh& m0,
    const DEBUG_PolygonMesh& m,
    const Tuple& h_tuple)
{
    std::array<long, 2> hids = {m0.id(h_tuple, PH), m0.id(m0.opp_halfedge(h_tuple), PH)};
    std::array<long, 2> hnids = {
        m0.id(m0.next_halfedge(h_tuple), PH),
        m0.id(m0.next_halfedge(m0.opp_halfedge(h_tuple)), PH)};
    std::array<long, 2> hpids = {
        m0.id(m0.prev_halfedge(h_tuple), PH),
        m0.id(m0.prev_halfedge(m0.opp_halfedge(h_tuple)), PH)};

    // Check distinct vertex loops
    if (m0.id(h_tuple, PV) != m0.id(m0.opp_halfedge(h_tuple), PV)) {
        for (long i = 0; i < 2; ++i) {
            std::vector<long> v0 = m0.get_vertex_id_loop(hids[i]);
            std::vector<long> v = m.get_vertex_id_loop_range(hpids[1 - i], hpids[i]);
            if ((v0.size() != v.size() + 1) || (!is_periodic_subsequence(v, v0, 1))) {
                return false;
            }
        }
    }
    // Check common vertex loops
    else {
        for (long i = 0; i < 2; ++i) {
            std::vector<long> v0 = m0.get_vertex_id_loop_range(hids[i], hids[1 - i]);
            std::vector<long> v = m.get_vertex_id_loop_range(hpids[1 - i], hpids[i]);
            if ((v0.size() != v.size() + 1) || (!is_periodic_subsequence(v, v0, 1))) {
                return false;
            }
        }
    }

    // Check distinct face loops
    if (m0.id(h_tuple, PF) != m0.id(m0.opp_halfedge(h_tuple), PF)) {
        for (long i = 0; i < 2; ++i) {
            std::vector<long> f0 = m0.get_face_id_loop(hids[i]);
            std::vector<long> f = m.get_face_id_loop(hnids[i]);
            if ((f0.size() != f.size() + 1) || (!is_periodic_subsequence(f, f0, 1))) {
                return false;
            }
        }
    }
    // Check common face loop
    else {
        std::vector<long> f0 = m0.get_face_id_loop(hids[0]);
        std::vector<long> f = m.get_face_id_loop(hnids[0]);
        if (f0.size() != f.size() + 2) {
            return false;
        }
        for (long i = 0; i < 2; ++i) {
            std::vector<long> f0 = m0.get_face_id_loop_range(hids[i], hids[1 - i]);
            std::vector<long> f = m.get_face_id_loop_range(hnids[i], hnids[1 - i]);
            if ((f0.size() != f.size() + 1) || (!is_periodic_subsequence(f, f0, 1))) {
                return false;
            }
        }
    }

    // Check the new faces are holes iff the original faces are holes
    if ((m.is_hole_face(m.tuple_from_id(PH, hnids[0])) != m0.is_hole_face(h_tuple)) ||
        (m.is_hole_face(m.tuple_from_id(PH, hnids[1])) !=
         m0.is_hole_face(m0.opp_halfedge(h_tuple)))) {
        return false;
    }

    return true;
}

TEST_CASE("join_vertex_operation", "[operations][join_vertex][polygon]")
{
    using namespace operations;

    DEBUG_PolygonMesh m0, m;
    OperationSettings<polygon_mesh::JoinVertex> op_settings;
    bool success;

    SECTION("join digon vertex")
    {
        m0 = digon();
        m = digon();
        Tuple h = m.tuple_from_id(PH, 0);
        polygon_mesh::JoinVertex op(m, h, op_settings);
        CHECK(op.precondition());
        success = op();
        Tuple v = op.return_tuple();

        CHECK(join_vertex_postcondition(m0, m, h));
        CHECK(m.get_all(PrimitiveType::Vertex).size() == 1);
        CHECK(m.get_all(PrimitiveType::Edge).size() == 1);
        CHECK(m.get_all(PrimitiveType::Face).size() == 2);
        CHECK(m.count_hole_faces() == 0);
        CHECK(m.count_interior_faces() == 2);
    }
    SECTION("join triangle vertex")
    {
        m0 = triangle();
        m = triangle();
        Tuple h = m.tuple_from_id(PH, 1);
        polygon_mesh::JoinVertex op(m, h, op_settings);
        CHECK(op.precondition());
        success = op();
        Tuple v = op.return_tuple();

        CHECK(join_vertex_postcondition(m0, m, h));
        CHECK(m.get_all(PrimitiveType::Vertex).size() == 2);
        CHECK(m.get_all(PrimitiveType::Edge).size() == 2);
        CHECK(m.get_all(PrimitiveType::Face).size() == 2);
        CHECK(m.count_hole_faces() == 0);
        CHECK(m.count_interior_faces() == 2);
    }
    SECTION("join vertex with common face")
    {
        m0 = torus();
        m = torus();
        Tuple h = m.tuple_from_id(PH, 0);
        polygon_mesh::JoinVertex op(m, h, op_settings);
        CHECK(op.precondition());
        success = op();
        Tuple v = op.return_tuple();

        CHECK(join_vertex_postcondition(m0, m, h));
        CHECK(m.get_all(PrimitiveType::Vertex).size() == 2);
        CHECK(m.get_all(PrimitiveType::Edge).size() == 1);
        CHECK(m.get_all(PrimitiveType::Face).size() == 1);
        CHECK(m.count_hole_faces() == 0);
        CHECK(m.count_interior_faces() == 1);
    }
    SECTION("join vertex with hole face")
    {
        m0 = triangle();
        polygon_mesh::MakeHole(m0, m0.tuple_from_id(PH, 0), {})();
        m = triangle();
        polygon_mesh::MakeHole(m, m.tuple_from_id(PH, 0), {})();

        Tuple h = m.tuple_from_id(PH, 0);
        polygon_mesh::JoinVertex op(m, h, op_settings);
        CHECK(op.precondition());
        success = op();
        Tuple v = op.return_tuple();

        CHECK(join_vertex_postcondition(m0, m, h));
        CHECK(m.get_all(PrimitiveType::Vertex).size() == 2);
        CHECK(m.get_all(PrimitiveType::Edge).size() == 2);
        CHECK(m.get_all(PrimitiveType::Face).size() == 2);
        CHECK(m.count_hole_faces() == 1);
        CHECK(m.count_interior_faces() == 1);
    }
    SECTION("join vertex with two hole faces")
    {
        m0 = triangle();
        polygon_mesh::MakeHole(m0, m0.tuple_from_id(PH, 0), {})();
        polygon_mesh::MakeHole(m0, m0.tuple_from_id(PH, 1), {})();
        m = triangle();
        polygon_mesh::MakeHole(m, m.tuple_from_id(PH, 0), {})();
        polygon_mesh::MakeHole(m, m.tuple_from_id(PH, 1), {})();

        Tuple h = m.tuple_from_id(PH, 0);
        polygon_mesh::JoinVertex op(m, h, op_settings);
        CHECK(op.precondition());
        success = op();
        Tuple v = op.return_tuple();

        CHECK(join_vertex_postcondition(m0, m, h));
        CHECK(m.get_all(PrimitiveType::Vertex).size() == 2);
        CHECK(m.get_all(PrimitiveType::Edge).size() == 2);
        CHECK(m.get_all(PrimitiveType::Face).size() == 2);
        CHECK(m.count_hole_faces() == 2);
        CHECK(m.count_interior_faces() == 0);
    }
    SECTION("join boundary and interior vertex")
    {
        m0 = grid();
        m = grid();

        Tuple h = m.halfedge_tuple_from_vertex_in_face(3, 2);
        polygon_mesh::JoinVertex op(m, h, op_settings);
        CHECK(op.precondition());
        success = op();
        Tuple v = op.return_tuple();

        CHECK(join_vertex_postcondition(m0, m, h));
        CHECK(m.get_all(PrimitiveType::Vertex).size() == 8);
        CHECK(m.get_all(PrimitiveType::Edge).size() == 11);
        CHECK(m.get_all(PrimitiveType::Face).size() == 5);
        CHECK(m.count_hole_faces() == 1);
        CHECK(m.count_interior_faces() == 4);
    }
    SECTION("invalid boundary vertex join")
    {
        m = annulus();

        Tuple h = m.next_halfedge(m.halfedge_tuple_from_vertex_in_face(2, 0));
        polygon_mesh::JoinVertex op(m, h, op_settings);
        CHECK(!op.precondition());
    }
    SECTION("invalid dual bubble join")
    {
        m = dual_bubble();

        Tuple h = m.tuple_from_id(PH, 0);
        polygon_mesh::JoinVertex op(m, h, op_settings);
        CHECK(!op.precondition());
    }
    SECTION("invalid bubble join")
    {
        m = monogon_with_leaf();

        Tuple h = m.tuple_from_id(PH, 0);
        polygon_mesh::JoinVertex op(m, h, op_settings);
        CHECK(!op.precondition());
    }
    SECTION("invalid bubble join")
    {
        m = bubble();

        Tuple h = m.tuple_from_id(PH, 0);
        polygon_mesh::JoinVertex op(m, h, op_settings);
        CHECK(!op.precondition());
    }

    CHECK(success);
    REQUIRE(m.is_connectivity_valid());
}
