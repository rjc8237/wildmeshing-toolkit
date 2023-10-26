set(SRC_FILES
	AtomicOperation.hpp
	AtomicOperation.cpp
	DeleteBubble.cpp
	DeleteBubble.hpp
	MakeEdge.cpp
	MakeEdge.hpp
	PolygonMeshOperation.cpp
	PolygonMeshOperation.hpp
	Splice.cpp
	Splice.hpp
)
target_sources(wildmeshing_toolkit PRIVATE ${SRC_FILES})
