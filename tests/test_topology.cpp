#include <wmtk/utils/tetmesh_topology_initialization.h>
#include <wmtk/utils/trimesh_topology_initialization.h>

#include <wmtk/Mesh.hpp>

#include <catch2/catch_test_macros.hpp>

#include <igl/readMESH.h>
#include <igl/readMSH.h>
#include <igl/read_triangle_mesh.h>

#include <stdlib.h>
#include <iostream>

#include <wmtk/utils/Logger.hpp>

using namespace wmtk;

namespace {
// Implicit opposite map for halfedges paired as [e] = {2*e, 2*e + 1}
long opp(long h)
{
    return ((h % 2) == 0) ? (h + 1) : (h - 1);
}
} // namespace

void check_polygon_mesh_edges(Eigen::Ref<const VectorXl> next, Eigen::Ref<const VectorXl> prev)
{
    REQUIRE(next.size() == prev.size());
    long n_halfedges = next.size();

    // Check that prev is a right and left inverse for next
    for (long hi = 0; hi < n_halfedges; ++hi) {
        CHECK(next[prev[hi]] == hi);
        CHECK(prev[next[hi]] == hi);
    }
}

void check_polygon_mesh_vertices(
    Eigen::Ref<const VectorXl> prev,
    Eigen::Ref<const VectorXl> to,
    Eigen::Ref<const VectorXl> out)
{
    REQUIRE(prev.size() == to.size());
    long n_halfedges = to.size();
    long n_vertices = out.size();

    // Check that to is invariant under circulation h -> prev[opp[h]]
    for (long hi = 0; hi < n_halfedges; ++hi) {
        CHECK(to[prev[opp(hi)]] == to[hi]);
    }

    // Check that out is a right inverse for h -> from[h] = to[opp[h]]
    for (long vi = 0; vi < n_vertices; ++vi) {
        CHECK(to[opp(out[vi])] == vi);
    }
}

void check_polygon_mesh_faces(
    Eigen::Ref<const VectorXl> next,
    Eigen::Ref<const VectorXl> he2f,
    Eigen::Ref<const VectorXl> f2he)
{
    REQUIRE(next.size() == he2f.size());
    long n_halfedges = he2f.size();
    long n_faces = f2he.size();

    // Check that he2f is invariant under next
    for (long hi = 0; hi < n_halfedges; ++hi) {
        CHECK(he2f[next[hi]] == he2f[hi]);
    }

    // Check that f2he is a right inverse for he2f
    for (long fi = 0; fi < n_faces; ++fi) {
        CHECK(he2f[f2he[fi]] == fi);
    }
}

TEST_CASE("topology_of_single_triangle", "[topology][2D]")
{
    Eigen::Matrix<long, 1, 3> F;
    F << 0, 1, 2;
    auto [FE, FF, VF, EF] = trimesh_topology_initialization(F);

    // std::cout << "F:\n" << F << std::endl;
    // std::cout << "FE:\n" << FE << std::endl;
    // std::cout << "FF:\n" << FF << std::endl;
    // std::cout << "VF:\n" << VF << std::endl;
    // std::cout << "EF:\n" << EF << std::endl;

    // 1. Test relationship between EF and FE
    for (int i = 0; i < EF.size(); ++i) {
        CHECK((FE.row(EF(i)).array() == i).any());
    }

    // 2. Test relationship between VF and F
    for (int i = 0; i < VF.size(); ++i) {
        CHECK((F.row(VF(i)).array() == i).any());
    }

    // 3. Test relationship between FF and FE
    for (int i = 0; i < FF.rows(); ++i) {
        for (int j = 0; j < 3; ++j) {
            long nb = FF(i, j);
            if (nb < 0) continue;

            CHECK((FF.row(nb).array() == i).any());

            if ((FF.row(nb).array() == i).any()) {
                int cnt = (FF.row(nb).array() == i).count();
                CHECK(cnt == 1);

                auto is_nb = (FF.row(nb).array() == i);
                for (int k = 0; k < 3; ++k) {
                    if (is_nb(k)) {
                        CHECK(FE(i, j) == FE(nb, k));
                    }
                }
            }
        }
    }
}

TEST_CASE("topology_of_two_triangles", "[topology][2D]")
{
    Eigen::Matrix<long, 2, 3> F;
    F << 0, 1, 2, 1, 3, 2;

    auto [FE, FF, VF, EF] = trimesh_topology_initialization(F);

    // std::cout << "F:\n" << F << std::endl;
    // std::cout << "FE:\n" << FE << std::endl;
    // std::cout << "FF:\n" << FF << std::endl;
    // std::cout << "VF:\n" << VF << std::endl;
    // std::cout << "EF:\n" << EF << std::endl;

    // 1. Test relationship between EF and FE
    for (int i = 0; i < EF.size(); ++i) {
        CHECK((FE.row(EF(i)).array() == i).any());
    }

    // 2. Test relationship between VF and F
    for (int i = 0; i < VF.size(); ++i) {
        CHECK((F.row(VF(i)).array() == i).any());
    }

    // 3. Test relationship between FF and FE
    for (int i = 0; i < FF.rows(); ++i) {
        for (int j = 0; j < 3; ++j) {
            long nb = FF(i, j);
            if (nb < 0) continue;

            CHECK((FF.row(nb).array() == i).any());

            if ((FF.row(nb).array() == i).any()) {
                int cnt = (FF.row(nb).array() == i).count();
                CHECK(cnt == 1);

                auto is_nb = (FF.row(nb).array() == i);
                for (int k = 0; k < 3; ++k) {
                    if (is_nb(k)) {
                        CHECK(FE(i, j) == FE(nb, k));
                    }
                }
            }
        }
    }
}

TEST_CASE("topology_of_complex_meshes", "[topology][2D]")
{
    Eigen::MatrixXd V;
    Eigen::Matrix<long, -1, -1> F;

    std::vector<std::string> names = {
        "/Octocat.obj",
        "/armadillo.obj",
        "/blub.obj",
        "/bunny.obj",
        "/circle.obj",
        "/fan.obj",
        "/sphere.obj",
        "/test_triwild.obj",
        "/hemisphere.obj"};

    for (auto name : names) {
        std::string path;
        path.append(WMTK_DATA_DIR);
        path.append(name);
        igl::read_triangle_mesh(path, V, F);

        auto [FE, FF, VF, EF] = trimesh_topology_initialization(F);

        // std::cout << "F:\n" << F << std::endl;
        // std::cout << "FE:\n" << FE << std::endl;
        // std::cout << "FF:\n" << FF << std::endl;
        // std::cout << "VF:\n" << VF << std::endl;
        // std::cout << "EF:\n" << EF << std::endl;

        // 1. Test relationship between EF and FE
        for (int i = 0; i < EF.size(); ++i) {
            CHECK((FE.row(EF(i)).array() == i).any());
        }

        // 2. Test relationship between VF and F
        for (int i = 0; i < VF.size(); ++i) {
            if (VF(i) < 0) continue;
            CHECK((F.row(VF(i)).array() == i).any());
        }

        // 3. Test relationship between FF and FE
        for (int i = 0; i < FF.rows(); ++i) {
            for (int j = 0; j < 3; ++j) {
                long nb = FF(i, j);
                if (nb < 0) continue;

                CHECK((FF.row(nb).array() == i).any());

                if ((FF.row(nb).array() == i).any()) {
                    int cnt = (FF.row(nb).array() == i).count();
                    CHECK(cnt == 1);

                    auto is_nb = (FF.row(nb).array() == i);
                    for (int k = 0; k < 3; ++k) {
                        if (is_nb(k)) {
                            // wmtk::logger().info("{} {} {} {}", i, j, nb, k);
                            CHECK(FE(i, j) == FE(nb, k));
                        }
                    }
                }
            }
        }
    }
}

TEST_CASE("topology_of_bubble", "[topology][polygon]")
{
    VectorXl next(2);
    next << 1, 0;
    std::vector<long> bnd_loops = {};
    auto [prev, to, he2f, f2he, out] = polygon_mesh_topology_initialization(next, bnd_loops);

    check_polygon_mesh_edges(next, prev);
    check_polygon_mesh_vertices(prev, to, out);
    check_polygon_mesh_faces(next, he2f, f2he);
}

TEST_CASE("topology_of_dual_bubble", "[topology][polygon]")
{
    VectorXl next(2);
    next << 0, 1;
    std::vector<long> bnd_loops = {};
    auto [prev, to, he2f, f2he, out] = polygon_mesh_topology_initialization(next, bnd_loops);

    check_polygon_mesh_edges(next, prev);
    check_polygon_mesh_vertices(prev, to, out);
    check_polygon_mesh_faces(next, he2f, f2he);
}

TEST_CASE("topology_of_single_triangle_polygon", "[topology][polygon]")
{
    Eigen::Matrix<long, 1, 3> F;
    F << 0, 1, 2;
    auto [next, prev, to, he2f, f2he, out, bnd_loops] = polygon_mesh_topology_initialization(F);

    check_polygon_mesh_edges(next, prev);
    check_polygon_mesh_vertices(prev, to, out);
    check_polygon_mesh_faces(next, he2f, f2he);
    CHECK(bnd_loops.size() == 1);
}

TEST_CASE("topology_of_two_triangle_polygons", "[topology][polygon]")
{
    Eigen::Matrix<long, 2, 3> F;
    F << 0, 1, 2, 1, 3, 2;
    auto [next, prev, to, he2f, f2he, out, bnd_loops] = polygon_mesh_topology_initialization(F);

    check_polygon_mesh_edges(next, prev);
    check_polygon_mesh_vertices(prev, to, out);
    check_polygon_mesh_faces(next, he2f, f2he);
    CHECK(bnd_loops.size() == 0);
}

TEST_CASE("topology_of_complex_polygon_meshes", "[topology][polygon]")
{
    Eigen::MatrixXd V;
    Eigen::Matrix<long, -1, -1> F;

    std::vector<std::string> names = {
        "/Octocat.obj",
        "/armadillo.obj",
        "/blub.obj",
        "/bunny.obj",
        "/circle.obj",
        "/fan.obj",
        "/sphere.obj",
        "/test_triwild.obj",
        "/hemisphere.obj"};

    for (auto name : names) {
        std::string path;
        path.append(WMTK_DATA_DIR);
        path.append(name);
        igl::read_triangle_mesh(path, V, F);
        auto [next, prev, to, he2f, f2he, out, bnd_loops] = polygon_mesh_topology_initialization(F);

        check_polygon_mesh_edges(next, prev);
        check_polygon_mesh_vertices(prev, to, out);
        check_polygon_mesh_faces(next, he2f, f2he);
    }
}


TEST_CASE("topology_of_two_adjacent_tets", "[topology][3D]")
{
    // Two tetrahedra are sharing one face
    // there are 7 unique faces and 9 unique edges

    Eigen::Matrix<long, 2, 4> T;
    T << 0, 1, 2, 3, 1, 2, 3, 4;

    auto [TE, TF, TT, VT, ET, FT] = tetmesh_topology_initialization(T);

    // 1. Test the maximum in TE and TF
    CHECK(TE.maxCoeff() == 9 - 1);
    CHECK(TF.maxCoeff() == 7 - 1);
    CHECK(TT.maxCoeff() == 1);

    // 2. Test the relationship between ET and TE
    for (int i = 0; i < ET.size(); ++i) {
        CHECK((TE.row(ET(i)).array() == i).any());
    }

    // 3. Test the relationship between FT and TF
    for (int i = 0; i < FT.size(); ++i) {
        CHECK((TF.row(FT(i)).array() == i).any());
    }

    // 4. Test the relationship between VT and T
    for (int i = 0; i < VT.size(); ++i) {
        if (VT(i) < 0) continue;
        CHECK((T.row(VT(i)).array() == i).any());
    }

    // 5. Test the relationship between TT and TF and TE
    for (int i = 0; i < TT.rows(); ++i) {
        for (int j = 0; j < 4; ++j) {
            long nb = TT(i, j);
            if (nb < 0) continue;

            CHECK((TT.row(nb).array() == i).any());

            if ((TT.row(nb).array() == i).any()) {
                int cnt = (TT.row(nb).array() == i).count();
                CHECK(cnt == 1);

                auto is_nb = (TT.row(nb).array() == i);
                for (int k = 0; k < 4; ++k) {
                    if (is_nb(k)) {
                        // wmtk::logger().info("{} {} {} {}", i, j, nb, k);
                        CHECK(TF(i, j) == TF(nb, k));
                    }
                }
            }
        }
    }
}

TEST_CASE("topology_of_two_independent_tets", "[topology][3D]")
{
    // Two tetrahedra not sharing anything
    // there are 8 unique faces and 12 unique edges

    Eigen::Matrix<long, 2, 4> T;
    T << 0, 1, 2, 3, 4, 5, 6, 7;

    auto [TE, TF, TT, VT, ET, FT] = tetmesh_topology_initialization(T);

    // 1. Test the maximum in TE and TF
    CHECK(TE.maxCoeff() == 12 - 1);
    CHECK(TF.maxCoeff() == 8 - 1);
    CHECK(TT.maxCoeff() == -1);

    // 2. Test the relationship between ET and TE
    for (int i = 0; i < ET.size(); ++i) {
        CHECK((TE.row(ET(i)).array() == i).any());
    }

    // 3. Test the relationship between FT and TF
    for (int i = 0; i < FT.size(); ++i) {
        CHECK((TF.row(FT(i)).array() == i).any());
    }

    // 4. Test the relationship between VT and T
    for (int i = 0; i < VT.size(); ++i) {
        if (VT(i) < 0) continue;
        CHECK((T.row(VT(i)).array() == i).any());
    }

    // 5. Test the relationship between TT and TF and TE
    for (int i = 0; i < TT.rows(); ++i) {
        for (int j = 0; j < 4; ++j) {
            long nb = TT(i, j);
            if (nb < 0) continue;

            CHECK((TT.row(nb).array() == i).any());

            if ((TT.row(nb).array() == i).any()) {
                int cnt = (TT.row(nb).array() == i).count();
                CHECK(cnt == 1);

                auto is_nb = (TT.row(nb).array() == i);
                for (int k = 0; k < 4; ++k) {
                    if (is_nb(k)) {
                        // wmtk::logger().info("{} {} {} {}", i, j, nb, k);
                        CHECK(TF(i, j) == TF(nb, k));
                    }
                }
            }
        }
    }
}

TEST_CASE("topology_of_tet_bunny", "[topology][3D]")
{
    Eigen::MatrixXd V;
    Eigen::Matrix<long, -1, -1> T, F;
    igl::readMESH(WMTK_DATA_DIR "/bunny.mesh", V, T, F);

    auto [TE, TF, TT, VT, ET, FT] = tetmesh_topology_initialization(T);

    // 1. Test the maximum in TE and TF
    CHECK(TE.maxCoeff() == (ET.size() - 1));
    CHECK(TF.maxCoeff() == (FT.size() - 1));

    // 2. Test the relationship between ET and TE
    for (int i = 0; i < ET.size(); ++i) {
        CHECK((TE.row(ET(i)).array() == i).any());
    }

    // 3. Test the relationship between FT and TF
    for (int i = 0; i < FT.size(); ++i) {
        CHECK((TF.row(FT(i)).array() == i).any());
    }

    // 4. Test the relationship between VT and T
    for (int i = 0; i < VT.size(); ++i) {
        if (VT(i) < 0) continue;
        CHECK((T.row(VT(i)).array() == i).any());
    }

    // 5. Test the relationship between TT and TF and TE
    for (int i = 0; i < TT.rows(); ++i) {
        for (int j = 0; j < 4; ++j) {
            long nb = TT(i, j);
            if (nb < 0) continue;

            CHECK((TT.row(nb).array() == i).any());

            if ((TT.row(nb).array() == i).any()) {
                int cnt = (TT.row(nb).array() == i).count();
                CHECK(cnt == 1);

                auto is_nb = (TT.row(nb).array() == i);
                for (int k = 0; k < 4; ++k) {
                    if (is_nb(k)) {
                        // wmtk::logger().info("{} {} {} {}", i, j, nb, k);
                        CHECK(TF(i, j) == TF(nb, k));
                    }
                }
            }
        }
    }
}
