set(SRC_FILES
    Invariant.hpp
    Invariant.cpp
    TriMeshInvariant.hpp
    TriMeshInvariant.cpp
    TetMeshInvariant.hpp
    TetMeshInvariant.cpp
    MeshInvariant.hpp
    MeshInvariant.cpp
    InvariantCollection.hpp
    InvariantCollection.cpp
    InteriorVertexInvariant.hpp
    InteriorVertexInvariant.cpp
    InteriorEdgeInvariant.hpp
    InteriorEdgeInvariant.cpp
    ValidTupleInvariant.hpp
    ValidTupleInvariant.cpp
    TriMeshLinkConditionInvariant.hpp
    TriMeshLinkConditionInvariant.cpp
    MultiMeshLinkConditionInvariant.hpp
    MultiMeshLinkConditionInvariant.cpp
    MaxEdgeLengthInvariant.hpp
    MaxEdgeLengthInvariant.cpp
    MinEdgeLengthInvariant.hpp
    MinEdgeLengthInvariant.cpp
    MinIncidentValenceInvariant.hpp
    MinIncidentValenceInvariant.cpp
    TriangleInversionInvariant.hpp
    TriangleInversionInvariant.cpp
    TodoInvariant.hpp
    TodoInvariant.cpp
)
target_sources(wildmeshing_toolkit PRIVATE ${SRC_FILES})
