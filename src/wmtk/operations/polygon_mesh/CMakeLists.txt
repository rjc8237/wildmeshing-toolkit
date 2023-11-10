set(SRC_FILES
	AtomicOperation.hpp
	AtomicOperation.cpp
	DeleteBubble.cpp
	DeleteBubble.hpp
	FillHole.cpp
	FillHole.hpp
	MakeEdge.cpp
	MakeEdge.hpp
	MakeHole.cpp
	MakeHole.hpp
	PolygonMeshOperation.cpp
	PolygonMeshOperation.hpp
	Splice.cpp
	Splice.hpp
)
target_sources(wildmeshing_toolkit PRIVATE ${SRC_FILES})
