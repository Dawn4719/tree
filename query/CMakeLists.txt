add_executable(IDE
        main.cpp
        deltagraph.h
        deltagraph.cpp
        ViLa.h
        ViLa.cpp
        ticktock.h
        staticcore.hpp
        LLAMA.hpp
        EvaluateQuery.h
        EvaluateQuery.cpp
        FilterVertices.cpp
        FilterVertices.h
        GenerateFilteringPlan.cpp
        GenerateFilteringPlan.h
        GenerateQueryPlan.cpp
        GenerateQueryPlan.h
        matchingcommand.cpp
        matchingcommand.h
        CECI.hpp
        IIT.hpp
        IITR.hpp
        graphDelta.cpp
        graphDelta.h
        type.h
)

find_library(ROCKSDB_LIB rocksdb)
find_package(Threads)
find_path(ROCKSDB_INCLUDE_DIR rocksdb/db.h
        PATHS "")
include_directories(${ROCKSDB_INCLUDE_DIR})

target_link_libraries(IDE Threads::Threads ${ROCKSDB_LIB} utility graph roaring)

