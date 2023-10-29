
#include <spdlog/spdlog.h>
#include <stdlib.h>
#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <iostream>
#include <tuple>
#include <wmtk/Tuple.hpp>
#include <wmtk/autogen/is_ccw.hpp>
#include <wmtk/autogen/local_switch_tuple.hpp>
#include <wmtk/autogen/tet_mesh/autogenerated_tables.hpp>
#include <wmtk/autogen/tet_mesh/is_ccw.hpp>
#include <wmtk/autogen/tet_mesh/local_id_table_offset.hpp>
#include <wmtk/autogen/tet_mesh/local_switch_tuple.hpp>
#include <wmtk/autogen/tri_mesh/autogenerated_tables.hpp>
#include <wmtk/autogen/tri_mesh/is_ccw.hpp>
#include <wmtk/autogen/tri_mesh/local_id_table_offset.hpp>
#include <wmtk/autogen/tri_mesh/local_switch_tuple.hpp>
#include <wmtk/multimesh/utils/find_local_switch_sequence.hpp>
#include <wmtk/multimesh/utils/local_switch_tuple.hpp>
#include <wmtk/multimesh/utils/transport_tuple.hpp>
#include <wmtk/utils/TupleInspector.hpp>
#include "tools/all_valid_local_tuples.hpp"

using namespace wmtk;
using namespace wmtk::tests;
TEST_CASE("tuple_autogen_find_all_local_switches", "[tuple]")
{
    for (wmtk::PrimitiveType pt :
         {PrimitiveType::Edge, PrimitiveType::Face, PrimitiveType::Tetrahedron}) {
        auto all = all_valid_local_tuples(pt);
        for (const Tuple& a : all) {
            for (const Tuple& b : all) {
                auto seq = wmtk::multimesh::utils::find_local_switch_sequence(a, b, pt);
                const Tuple myb = wmtk::multimesh::utils::local_switch_tuples(pt, a, seq);

                spdlog::warn(
                    "{}=>{} == {}",
                    wmtk::utils::TupleInspector::as_string(a),
                    wmtk::utils::TupleInspector::as_string(myb),
                    wmtk::utils::TupleInspector::as_string(b));
                CHECK(myb == b);
            }
        }
    }
}
