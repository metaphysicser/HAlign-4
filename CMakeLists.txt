cmake_minimum_required(VERSION 3.18)
project(HAlign4)

set(CMAKE_CXX_STANDARD 20)
option(HALIGN4_NATIVE_ARCH "Enable -march=native" OFF)

# 如果用户没有指定 CMAKE_BUILD_TYPE，则默认设置为 Release，以便默认开启优化构建
if(NOT DEFINED CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build (default: Release)" FORCE)
    message(STATUS "CMAKE_BUILD_TYPE not set -> defaulting to: ${CMAKE_BUILD_TYPE}")
else()
    message(STATUS "CMAKE_BUILD_TYPE is: ${CMAKE_BUILD_TYPE}")
endif()


include_directories(
        ${PROJECT_SOURCE_DIR}/include
        ${PROJECT_SOURCE_DIR}/extern/CLI11/include
        ${PROJECT_SOURCE_DIR}/extern/spdlog/include
        ${PROJECT_SOURCE_DIR}/extern/cereal/include
        ${PROJECT_SOURCE_DIR}/extern/WFA2-lib
)

add_executable(halign4
        src/halign4.cpp
        src/utils/preprocess.cpp
        src/utils/seq_io.cpp
        src/utils/cmd.cpp
        src/utils/file_io.cpp
        src/mash/hash.cpp
        src/mash/xxhash.c
        src/mash/mash.cpp
        src/consensus/selector.cpp
        src/consensus/consensus.cpp
        src/align/ref_aligner.cpp
        src/align/ksw2_extz2_sse.c
        src/align/psw_gg3_sse.c
        src/align/cigar.cpp
        src/align/align.cpp
        src/seed/minimizer.cpp
        src/seed/anchor.cpp

)

find_package(CURL)
if (CURL_FOUND)
    target_link_libraries(halign4 PRIVATE CURL::libcurl)
    target_compile_definitions(halign4 PRIVATE HALIGN4_HAVE_LIBCURL=1)
else()
    message(STATUS "libcurl not found; downloadFile will fall back to calling curl/wget at runtime")
endif()


find_package(ZLIB)
if (ZLIB_FOUND)
    target_link_libraries(halign4 PRIVATE ZLIB::ZLIB)
    target_compile_definitions(halign4 PRIVATE HALIGN4_HAVE_ZLIB=1)
else()
    message(STATUS "zlib not found; compression support will be disabled or limited")
endif()

find_package(OpenMP REQUIRED)

# ---------------------------------------------------------------------------
# WFA2-lib（Wavefront Alignment）子项目裁剪：只编译核心库，不编译测试/benchmark
# ---------------------------------------------------------------------------
set(BUILD_TESTING OFF CACHE BOOL "Disable tests for subprojects" FORCE)
set(CMAKE_DISABLE_FIND_PACKAGE_PkgConfig ON CACHE BOOL "Avoid PkgConfig dependency for WFA2-lib" FORCE)

# WFA2-lib 的 OPENMP 选项默认 OFF。但我们的主工程已启用 OpenMP，
# 为了让 WFA2-lib 正确打开 WFA_PARALLEL，这里显式打开。
set(OPENMP ON CACHE BOOL "Enable OpenMP for WFA2-lib" FORCE)

# 关键：只编译核心库，关闭 WFA2-lib 自带 benchmark/tests。
# （这两个开关由我们在 extern/WFA2-lib/CMakeLists.txt 中补充，默认 ON；这里只在主工程里强制关闭。）
set(WFA2LIB_BUILD_BENCHMARK OFF CACHE BOOL "Do not build WFA2-lib benchmark" FORCE)
set(WFA2LIB_BUILD_TESTS OFF CACHE BOOL "Do not register WFA2-lib tests" FORCE)

add_subdirectory(extern/WFA2-lib)

target_link_libraries(halign4 PRIVATE OpenMP::OpenMP_CXX wfa2::wfa2cpp_static)

# 在 Release 时添加激进的优化标志：GCC/Clang 与 MSVC 分别处理
if (MSVC)
    # MSVC 优化：/O2 开启优化，/arch:AVX2 请求 AVX2 指令集（可根据需要改为 AVX 或 AVX512），/GL 全局优化
    target_compile_options(halign4 PRIVATE
        $<$<CONFIG:Release>:/O2>
        $<$<CONFIG:Release>:/arch:AVX2>
        $<$<CONFIG:Release>:/GL>
    )
    # 链接时启用链接时期优化（LTCG）
    target_link_options(halign4 PRIVATE
        $<$<CONFIG:Release>:/LTCG>
    )
else()
    # GCC / Clang: -O3, -march=native, -flto, -pipe
    target_compile_options(halign4 PRIVATE
        $<$<CONFIG:Release>:-O3>
            $<$<AND:$<CONFIG:Release>,$<BOOL:${HALIGN4_NATIVE_ARCH}>>:-march=native>

    )
endif()


target_compile_definitions(halign4 PRIVATE
        $<$<CONFIG:Debug>:_DEBUG>
)
