# CMakeList.txt : CMake project for interface_simulink, include source and define
# project specific logic here.
#
cmake_minimum_required (VERSION 3.8)

project ("interface_simulink")

add_executable (lib_call "main_lib_call.cpp" "rt_defines.h" "rt_nonfinite.h" "rtGetInf.h" "rtGetNaN.h" "uavNavigationRL2D_2022a_icra.h" "uavNavigationRL2D_2022a_icra_private.h" "uavNavigationRL2D_2022a_icra_types.h")

target_link_libraries(lib_call -ldl) #https://stackoverflow.com/questions/11783932/how-do-i-add-a-linker-or-compile-flag-in-a-cmake-file

# Add library
add_compile_definitions(RT_MALLOC MODEL=uavNavigationRL2D_2022a_icra NUMST=2 NCSTATES=14 HAVESTDIO RT USE_RTMODEL CLASSIC_INTERFACE=0 ALLOCATIONFCN=1 TID01EQ=1 MAT_FILE=0 ONESTEPFCN=1 TERMFCN=1 MULTI_INSTANCE_CODE=1 INTEGER_CODE=0 MT=0)

#SET(GCC_COVERAGE_LINK_FLAGS    "-z defs")

include(GNUInstallDirs)
add_library(interface_simulink SHARED "wrapper_routine_RL2D.cpp" "ext_work.h" "multiword_types.h" "rt_logging.h" "rtmodel.h" "rtw_continuous.h" "rtw_solver.h" "rtwtypes.h" "tmwtypes.h" "rt_defines.h" "rt_nonfinite.cpp" "rt_nonfinite.h" "rtGetInf.cpp" "rtGetInf.h" "rtGetNaN.cpp" "rtGetNaN.h" "uavNavigationRL2D_2022a_icra.cpp" "uavNavigationRL2D_2022a_icra.h" "uavNavigationRL2D_2022a_icra_data.cpp" "uavNavigationRL2D_2022a_icra_private.h" "uavNavigationRL2D_2022a_icra_types.h")
set_target_properties(interface_simulink PROPERTIES
    VERSION 1
    SOVERSION 1
    PUBLIC_HEADER wrapper_api.h)

configure_file(interface_simulink.pc.in interface_simulink.pc @ONLY)

target_include_directories(interface_simulink PRIVATE .)

install(TARGETS interface_simulink
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

install(FILES ${CMAKE_BINARY_DIR}/interface_simulink.pc
    DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/pkgconfig)

#SET(CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} ${GCC_COVERAGE_LINK_FLAGS}")