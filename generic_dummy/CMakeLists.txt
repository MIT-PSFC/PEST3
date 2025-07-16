# $Id$

include_directories(${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/include
)

set(LIB_FILES
trread.f90
trxplib.f90
)

add_library(generic_dummy ${LIB_FILES})

