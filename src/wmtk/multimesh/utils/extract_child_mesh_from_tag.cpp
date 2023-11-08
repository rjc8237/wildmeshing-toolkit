#include "extract_child_mesh_from_tag.hpp"
#include <map>
#include <wmtk/EdgeMesh.hpp>
#include <wmtk/Mesh.hpp>
#include <wmtk/MultiMeshManager.hpp>
#include <wmtk/PointMesh.hpp>
#include <wmtk/Primitive.hpp>
#include <wmtk/TetMesh.hpp>
#include <wmtk/TriMesh.hpp>


namespace wmtk::multimesh::utils {


void extract_and_register_child_mesh_from_tag(
    TriMesh& m,
    const std::string& tag,
    const long& tag_value,
    const PrimitiveType& pt)
{
    assert(m.top_simplex_type() >= pt);
    auto tag_handle = m.get_attribute_handle<long>(tag, pt);
    auto tags = m.create_const_accessor(tag_handle);

    std::vector<Tuple> tagged_tuples;
    for (auto t : m.get_all(pt)) {
        if (tags.const_scalar_attribute(t) == tag_value) {
            tagged_tuples.emplace_back(t);
        }
    }

    switch (pt) {
    case PrimitiveType::Vertex: throw("not implemented");
    case PrimitiveType::Edge: {
        std::map<long, long> parent_to_child_vertex_map;
        long child_vertex_count = 0;

        RowVectors2l edge_mesh_matrix;
        edge_mesh_matrix.resize(tagged_tuples.size(), 2);

        for (long i = 0; i < tagged_tuples.size(); ++i) {
            const long v1 = m.id(tagged_tuples[i], PrimitiveType::Vertex);
            const long v2 = m.id(m.switch_vertex(tagged_tuples[i]), PrimitiveType::Vertex);

            // check and add v1, v2 to the vertex map
            if (parent_to_child_vertex_map.find(v1) == parent_to_child_vertex_map.end()) {
                parent_to_child_vertex_map[v1] = child_vertex_count;
                ++child_vertex_count;
            }
            if (parent_to_child_vertex_map.find(v2) == parent_to_child_vertex_map.end()) {
                parent_to_child_vertex_map[v2] = child_vertex_count;
                ++child_vertex_count;
            }

            // add edge to matrix
            edge_mesh_matrix(i, 0) = parent_to_child_vertex_map[v1];
            edge_mesh_matrix(i, 1) = parent_to_child_vertex_map[v2];
        }

        std::shared_ptr<EdgeMesh> child_ptr = std::make_shared<EdgeMesh>();
        auto& child = *child_ptr;
        child.initialize(edge_mesh_matrix);

        std::vector<std::array<Tuple, 2>> child_to_parent_map(tagged_tuples.size());
        assert(tagged_tuples.size() == child.capacity(PrimitiveType::Edge));

        auto edgemesh_edge_tuples = child.get_all(PrimitiveType::Edge);

        for (long i = 0; i < tagged_tuples.size(); ++i) {
            child_to_parent_map[i] = {{edgemesh_edge_tuples[i], tagged_tuples[i]}};
        }

        m.register_child_mesh(child_ptr, child_to_parent_map);
        break;
    }
    case PrimitiveType::Face: throw("not implemented");
    case PrimitiveType::Tetrahedron: throw("cannot register tetmesh on trimesh");
    default: throw("invalid child mesh type");
    }
}


} // namespace wmtk::multimesh::utils