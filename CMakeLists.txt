# cmake settings for qAlgorithms

cmake_minimum_required(VERSION 3.12)
project(qAlgorithms VERSION 0.2 LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Let CMake's selected configuration control optimisation and debug symbols.
if(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Build type" FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Debug Release RelWithDebInfo MinSizeRel)
endif()

option(QALGORITHMS_NATIVE_ARCH "Optimise for the build machine's CPU (not portable)" OFF)
option(QALGORITHMS_USE_ZLIB_NG "Use the native zlib-ng API instead of standard zlib" OFF)
include(CheckCXXCompilerFlag)

# Optional diagnostics vary between compiler versions. Only add supported flags.
function(qalgo_add_supported_option flag)
    string(MAKE_C_IDENTIFIER "QALGORITHMS_HAS_${flag}" flag_variable)
    check_cxx_compiler_flag("${flag}" ${flag_variable})
    if(${flag_variable})
        add_compile_options("${flag}")
    endif()
endfunction()

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
foreach(flag IN ITEMS
    -fno-math-errno            # Disables math functions setting errno, might result in speedup
    -fno-signed-zeros          # no difference between +0 and -0, might increase performance
    -fstrict-flex-arrays=3     # Size of arrays at the end of a struct must be known at compile time
    -Wall                      # Enables all recommended warnings
    -Wextra                    # Enables additional warnings
    -Wpedantic                 # Warns if not standard compliant
    -Wuninitialized            # Warns about uninitialized variables
    -Wformat=2                 # Check calls to printf and scanf, etc., to make sure that the arguments have types appropriate to the format string
    -Wformat-truncation      # Warns about potential truncation of output in snprintf etc.
    -Wimplicit-fallthrough     # Warns on missing break in switch case block
    -Wshadow                   # Warns when variable with same name as already declared variable is declared in lower scope
    -Wnon-virtual-dtor         # Warns if a class has a non-virtual destructor
    -Wcast-align               # Warns about pointer casts that may misalign data
    -Wunused                   # Warns about unused variables or functions
    -Wnull-dereference         # Warns about dereferencing null pointers
    -Wstrict-overflow=4        # Warns about cases where the compiler optimizes based on the assumption that signed overflow does not occur
    -fdiagnostics-color=always # Enables colored output for diagnostics
    -Wmissing-braces           # list initialisation with incorrect braces
    -Wredundant-decls          # warns about redundant declarations of the same entity
    -Warray-bounds           # compiler checks array bounds where possible
    -Wunreachable-code         # unreachable branches in code
    -Wenum-compare             # only compare enums of different types
    -Wpointer-arith
    -Waddress-of-packed-member
)
    qalgo_add_supported_option("${flag}")
endforeach()
add_compile_definitions($<$<CONFIG:Debug>:_GLIBCXX_ASSERTIONS>)
add_compile_options($<$<CONFIG:Debug>:-fno-omit-frame-pointer>)
add_compile_options($<$<CONFIG:Debug>:-ftrapv>) # Trap signed integer overflow in debug builds.
if(QALGORITHMS_NATIVE_ARCH)
    qalgo_add_supported_option(-march=native)
endif()
elseif(MSVC)
    add_compile_options(/W4)
endif()

# narrowing float conversion warning is different between clang and gcc
message(STATUS "Compiling with ${CMAKE_CXX_COMPILER_ID}")
if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    qalgo_add_supported_option(-Wimplicit-float-conversion)
    message(STATUS "Clang used")
elseif (CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    foreach(flag IN ITEMS
        -Wfloat-conversion
        -Wvector-operation-performance
        -Wduplicated-cond          # Warns if both arguments of conditional statement are the same (int a = (bool) ? 1 : 1)
        -Wlogical-op               # Warns about suspicious logical operations like (a < 0 && a < 0)
        -Wuseless-cast             # Warns casting expression to its own type
        -Wchanges-meaning          # @todo what does this do exactly?
        )
        qalgo_add_supported_option("${flag}")
    endforeach()
elseif(MSVC)
    message(STATUS "MSVC does not support a warning for float -> double conversions")
else()
    message(STATUS "Unknown compiler (${CMAKE_CXX_COMPILER_ID}): not adding float -> double conversion warning")
endif()

# flags for static compile 
# set(BUILD_SHARED_LIBS OFF)
# set(CMAKE_EXE_LINKER_FLAGS "-static")

# warnings that produce too much noise when building - @todo add these in eventually
# add_compile_options(-Wdouble-promotion -Wconversion -Weffc++)
# add_compile_options(-Wconversion)

# add_compile_options(-fsanitize=address) # @todo find a way to use this (not supported under mingw)
# add_link_options($<$<CONFIG:Debug>:-fsanitize=address>)

# lets make the directories into variables so we dont have to type CMAKE_CURRENT_SOURCE_DIR everywhere... 
set(TESTS_DIR "${CMAKE_CURRENT_SOURCE_DIR}/test")        # Directory containing test files and shared test resources
set(EXTERNAL_DIR "${CMAKE_CURRENT_SOURCE_DIR}/external") # Directory containing external headers/libs
set(SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src")        # Directory containing the Source Files
set(INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/include")   # Directory containing the Include Files

# set some warnings only inside the actual qalgo part of the code due to noise from libraries
set_source_files_properties(${SOURCE_DIR} PROPERTIES COMPILE_FLAGS 
    "-Wfloat-equal -Wsign-conversion -Wduplicated-branches"          
    # float == float emits a warning
    # conversion between signed and unsigned types
    # Warn when an if-else has identical branches. This warning also warns for conditional operators. Included here since the simdutf library is misdiagnosed by this a lot
)
set_source_files_properties(${TESTS_DIR} PROPERTIES COMPILE_FLAGS 
    "-Wfloat-equal -Wsign-conversion -Wduplicated-branches")

# Set include directories
include_directories(include "${EXTERNAL_DIR}/simdutf/" "${EXTERNAL_DIR}/cephes" "${EXTERNAL_DIR}/libcerf")  # @todo this is terrible practice. Lets not do it

# Source files
# @todo we currently list external files manually. Thats not great. 
### IMPORTANT: While using g++ to compile, all source files need to end in ".cpp". If a file ends in ###
### ".c", it will not be compiled and cause linker errors. This does not generate an error message!  ###
file(GLOB_RECURSE SOURCES 
    "${SOURCE_DIR}/*.cpp" 
    # "${EXTERNAL_DIR}/CDFlib/cdflib.cpp"
    "${EXTERNAL_DIR}/simdutf/simdutf.cpp" 
    "${EXTERNAL_DIR}/cephes/cephes.cpp"
    "${EXTERNAL_DIR}/libcerf/libcerf_reduced.cpp"
    "${EXTERNAL_DIR}/fastexp.h" 
    "${EXTERNAL_DIR}/exptables.h")

list(FILTER SOURCES EXCLUDE REGEX "qalgorithms_main\\.cpp") # Remove test/main file(s) from the library list

# Now we make a library from our sources that we can then link main.cpp and all the tests against
add_library(qalgo_lib ${SOURCES})
# Expose include dirs to consumers (app + tests)
target_include_directories(qalgo_lib PUBLIC
    ${INCLUDE_DIR}
    ${TESTS_DIR}/shared
    ${EXTERNAL_DIR}
)


# PROGRAM SOURCE
add_executable(${PROJECT_NAME} "${SOURCE_DIR}/qalgorithms_main.cpp")

# Keep the selected header and linked library in agreement. Users can pass
# ZLIB_ROOT or CMAKE_PREFIX_PATH for a non-standard installation, including MSYS2.
if(QALGORITHMS_USE_ZLIB_NG)
    find_package(zlib-ng CONFIG REQUIRED)
    target_compile_definitions(qalgo_lib PUBLIC QALGORITHMS_USE_ZLIB_NG)
    target_link_libraries(qalgo_lib PUBLIC zlib-ng::zlib)
else()
    find_package(ZLIB REQUIRED)
    target_link_libraries(qalgo_lib PUBLIC ZLIB::ZLIB)
endif()
target_link_libraries(${PROJECT_NAME} PRIVATE qalgo_lib)

# TESTS
include(CTest)
if(BUILD_TESTING)

file(GLOB_RECURSE TESTFILES "${TESTS_DIR}/*.cpp")
list(FILTER TESTFILES EXCLUDE REGEX ".*/test/shared/.*") # lets just assume that no one is stupid enough to have a second test/shared folder in this project

# Make the shared test code into a small library
file(GLOB_RECURSE TEST_SHARED_SOURCES "${TESTS_DIR}/shared/*.cpp")
add_library(test_shared_lib STATIC ${TEST_SHARED_SOURCES})
# @todo not sure which dirs we need to include here. 
target_include_directories(test_shared_lib PUBLIC ${TESTS_DIR}/shared ${INCLUDE_DIR} ${EXTERNAL_DIR})

# Now lets make the tests
foreach(src IN LISTS TESTFILES)
    # create a unique target name from relative path under tests/
    file(RELATIVE_PATH rel "${TESTS_DIR}" "${src}")
    string(REPLACE "/" "_" rel_safe "${rel}") # replace / with _
    string(REGEX REPLACE "\\.cpp$" "" rel_safe "${rel_safe}") # remove .cpp extension from filename
    set(target_name "${rel_safe}")
    message(NEW_TEST: "${target_name}")

    add_executable(${target_name} "${src}")
    set_target_properties(${target_name} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/tests_bin)
    target_link_libraries(${target_name} PRIVATE qalgo_lib test_shared_lib)
    add_test(NAME ${target_name} COMMAND ${target_name})
endforeach()

add_test(NAME test_qpeaks_regressions COMMAND test_qpeaks --regressions-only)
add_test(NAME test_qpeaks_pt03 COMMAND test_qpeaks --legacy-pt03)
add_test(NAME test_qpeaks_pt04 COMMAND test_qpeaks --legacy-pt04)

find_package(Python3 COMPONENTS Interpreter QUIET)
if(Python3_Interpreter_FOUND)
    add_test(NAME test_mzml_cli
        COMMAND "${Python3_EXECUTABLE}" "${TESTS_DIR}/integration/test_mzml_cli.py"
                "$<TARGET_FILE:qAlgorithms>")
    set_tests_properties(test_mzml_cli PROPERTIES TIMEOUT 120)
    add_test(NAME test_folding_cli
        COMMAND "${Python3_EXECUTABLE}" "${TESTS_DIR}/integration/test_folding_cli.py"
                "$<TARGET_FILE:qAlgorithms>")
    set_tests_properties(test_folding_cli PROPERTIES TIMEOUT 120)
endif()
endif()


# Optional standalone research executables: no production qPeaks integration.
option(QALGORITHMS_ENABLE_HYPERBOLIC_RESEARCH
       "Build isolated hyperbolic peak-running regression tests and benchmark" OFF)
if(QALGORITHMS_ENABLE_HYPERBOLIC_RESEARCH)
    add_executable(test_hyperbolic_running_research
        experiments/hyperbolic_running/test.cpp)
    add_executable(bench_hyperbolic_running_research
        experiments/hyperbolic_running/bench.cpp)
    target_link_libraries(bench_hyperbolic_running_research PRIVATE qalgo_lib)
    add_executable(real_hyperbolic_qtof_research experiments/hyperbolic_running/real_qtof.cpp)
    target_link_libraries(real_hyperbolic_qtof_research PRIVATE qalgo_lib)
    if(BUILD_TESTING)
        add_test(NAME test_hyperbolic_running_research
                 COMMAND test_hyperbolic_running_research)
    endif()
endif()

if(BUILD_TESTING)
    target_compile_definitions(test_peak_model PRIVATE STAGE_A)
endif()

option(QALGORITHMS_ENABLE_PEAK_BINDING_BENCH "Build isolated model-binding benchmark" OFF)
if(QALGORITHMS_ENABLE_PEAK_BINDING_BENCH)
    add_executable(bench_peak_binding experiments/peak_model_binding/bench.cpp)
    target_link_libraries(bench_peak_binding PRIVATE qalgo_lib)
    add_executable(bench_peak_binding_extra experiments/peak_model_binding/bench.cpp)
    target_compile_definitions(bench_peak_binding_extra PRIVATE STAGE_B)
    target_link_libraries(bench_peak_binding_extra PRIVATE qalgo_lib)
endif()

option(QALGORITHMS_ENABLE_PEAK_INFERENCE_BENCH "Build selected-window inference benchmark" OFF)
if(QALGORITHMS_ENABLE_PEAK_INFERENCE_BENCH)
    add_executable(bench_peak_inference experiments/peak_model_inference/bench.cpp)
    target_link_libraries(bench_peak_inference PRIVATE qalgo_lib)
endif()

option(QALGORITHMS_ENABLE_PEAK_SEARCH_BENCH "Build declared-family search benchmark" OFF)
if(QALGORITHMS_ENABLE_PEAK_SEARCH_BENCH)
    add_executable(bench_peak_search experiments/peak_model_search/bench.cpp)
    target_link_libraries(bench_peak_search PRIVATE qalgo_lib)
endif()

option(QALGORITHMS_ENABLE_REAL_PROFILE_AUDIT "Build isolated real-profile residual audit" OFF)
if(QALGORITHMS_ENABLE_REAL_PROFILE_AUDIT)
    add_executable(real_peak_audit experiments/peak_model_audit/collect.cpp)
    target_link_libraries(real_peak_audit PRIVATE qalgo_lib)
    if(BUILD_TESTING AND Python3_Interpreter_FOUND)
        add_test(NAME test_real_profile_residual_reference
            COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/experiments/peak_model_audit/test_analyze.py")
    endif()
endif()
