set(UTILS_SRC
    CLI11.hpp
    globals.cpp
    globals.h
    pod.h
    types.h
    utils.h
)

add_library(utils SHARED
        ${UTILS_SRC})