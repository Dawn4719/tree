set(UTILITY_SRC
        graphoperations.cpp
        graphoperations.h
        commandparser.cpp
        commandparser.h
        computesetintersection.h
        computesetintersection.cpp
        bitsetoperation.cpp
        bitsetoperation.h
        QFliter.h
        han/intersection_algos.cpp
        han/utils/util.cpp
utils.h
CLI11.hpp
globals.h
globals.cpp
types.h)

add_library(utility SHARED
        ${UTILITY_SRC})