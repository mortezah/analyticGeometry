#!/bin/bash -ex
export CMAKE_PREFIX_PATH= $SCORECPATH/lib/cmake/SCOREC

SCOREC=${PUMI_INSTALL_DIR}
cmake \
-DCMAKE_C_COMPILER="mpicc" \
-DCMAKE_CXX_COMPILER=mpicxx \
-DCMAKE_C_FLAGS=" -g -O2 -Wall -Wextra " \
-DCMAKE_CXX_FLAGS=" -g -O2 -Wall -Wextra " \
-DSCOREC_PREFIX=${SCOREC}
#-DCMAKE_CXX_FLAGS="-std=c++11" \

