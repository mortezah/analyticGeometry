cmake_minimum_required(VERSION 3.0)

#find_package(SCOREC REQUIRED PATHS ${SCOREC_DIR} NO_DEFAULT_PATH)
find_package(SCOREC REQUIRED)
#find_package(MPI_REQUIRED)
add_executable(main main.cpp)
target_link_libraries(main SCOREC::core)
target_link_libraries(main libmpi.so)

