//
// Created by Yixin Hu on 11/3/21.
//

#include "TetMesh.h"

void wmtk::TetMesh::split_edge(const Tuple &loc0) {
    std::shared_ptr<InfoCache> info = std::make_shared<InfoCache>();
    if (!split_before(loc0, info))
        return;

    // backup of everything
    auto loc1 = loc0;
    int v1_id = loc1.get_vid();
    auto loc2 = loc1.switch_vertex(*this);
    int v2_id = loc2.get_vid();
    //
    auto n12_t_ids = set_intersection(m_vertex_connectivity[v1_id].m_conn_tets,
                                      m_vertex_connectivity[v2_id].m_conn_tets);
    std::vector<size_t> n12_v_ids;
    for (size_t t_id: n12_t_ids) {
        for (int j = 0; j < 4; j++)
            n12_v_ids.push_back(m_tet_connectivity[t_id][j]);
    }
    vector_unique(n12_v_ids);
    std::vector<std::pair<size_t, TetrahedronConnectivity>> old_tets(n12_t_ids.size());
    std::vector<std::pair<size_t, VertexConnectivity>> old_vertices(n12_v_ids.size());
    for (size_t t_id: n12_t_ids)
        old_tets.push_back(std::make_pair(t_id, m_tet_connectivity[t_id]));
    for (size_t v_id: n12_v_ids)
        old_vertices.push_back(std::make_pair(v_id, m_vertex_connectivity[v_id]));

    // update connectivity
    int v_id = find_next_empty_slot_v();
    std::vector<size_t> new_t_ids;
    for (size_t t_id: n12_t_ids) {
        size_t new_t_id = find_next_empty_slot_t();
        new_t_ids.push_back(new_t_id);
        //
        int j = m_tet_connectivity[t_id].find(v1_id);
        m_tet_connectivity[new_t_id] = m_tet_connectivity[t_id];
        m_tet_connectivity[new_t_id][j] = v_id;
        //
        m_vertex_connectivity[v_id].m_conn_tets.push_back(t_id);
        m_vertex_connectivity[v_id].m_conn_tets.push_back(new_t_id);
        //
        for (int j = 0; j < 4; j++) {
            if (m_tet_connectivity[t_id][j] != v1_id && m_tet_connectivity[t_id][j] != v2_id)
                m_vertex_connectivity[m_tet_connectivity[t_id][j]].m_conn_tets.push_back(new_t_id);
        }
        //
        j = m_tet_connectivity[t_id].find(v2_id);
        m_tet_connectivity[t_id][j] = v_id;
        //
        vector_erase(m_vertex_connectivity[v1_id].m_conn_tets, t_id);
    }

    // possibly call the resize_attributes
    resize_attributes(m_vertex_connectivity.size(), m_tet_connectivity.size() * 6,
                      m_tet_connectivity.size() * 4, m_tet_connectivity.size());

    Tuple loc(v_id, 0, 0, old_tets.front().first);
    if (!split_after(loc, info)) {
        m_vertex_connectivity[v_id].m_is_removed = true;
        for (int t_id: new_t_ids)
            m_tet_connectivity[t_id].m_is_removed = true;

        //todo copy back old tets/vertices

        return;
    }

    // call invariants on all entities
    if (false) // if any invariant fails
    {
        // undo changes
        return;
    }

}