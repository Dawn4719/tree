set(ROARING_SRC
    roaring.c
    roaring.h
    roaring.hh
)

add_library(roaring SHARED
        ${ROARING_SRC})