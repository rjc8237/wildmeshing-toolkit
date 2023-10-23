#pragma once

// DO NOT MODIFY, autogenerated from the /scripts directory

namespace wmtk::autogen::tet_mesh {
    // lvids for a specfic leid
    extern const long auto_3d_edges[6][2];

    // lvids for a specific lfid
    extern const long auto_3d_faces[4][3];

    // vertex switch's (tuple_index -> [lvid,leid,lfid])
    extern const long auto_3d_table_vertex[96][3];

    // edge switch's (tuple_index -> [lvid,leid,lfid])
    extern const long auto_3d_table_edge[96][3];

    // face switch's (tuple_index -> [lvid,leid,lfid])
    extern const long auto_3d_table_face[96][3];

    // (tuple_index) -> is_ccw
    extern const long auto_3d_table_ccw[96];

    // lvid -> a ccw [lvid,leid,lfid]
    extern const long auto_3d_table_complete_vertex[4][3];

    // evid -> a ccw [lvid,leid,lfid]
    extern const long auto_3d_table_complete_edge[6][3];

    // fvid -> a ccw [lvid,leid,lfid]
    extern const long auto_3d_table_complete_face[4][3];

    
}