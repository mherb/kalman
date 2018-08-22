#  Eigen3_INCLUDE_DIR
#  Eigen3_FOUND

find_path( Eigen3_INCLUDE_DIR Eigen/Core
    PATH_SUFFIXES eigen3
)

include(LibFindMacros)
set(Eigen3_PROCESS_INCLUDES Eigen3_INCLUDE_DIR)
libfind_process(Eigen3)
