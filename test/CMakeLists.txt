cmake_minimum_required(VERSION 3.28)
project(HAlign4Tests LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


include(CTest)  # 提供 ctest 基础设施与 add_test :contentReference[oaicite:1]{index=1}

# 如果用户没有指定 CMAKE_BUILD_TYPE，则默认设置为 Release，以便默认开启优化构建
if(NOT DEFINED CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build (default: Release)" FORCE)
    message(STATUS "CMAKE_BUILD_TYPE not set -> defaulting to: ${CMAKE_BUILD_TYPE}")
else()
    message(STATUS "CMAKE_BUILD_TYPE is: ${CMAKE_BUILD_TYPE}")
endif()
get_filename_component(HALIGN4_ROOT "${CMAKE_CURRENT_LIST_DIR}/.." ABSOLUTE)


find_package(OpenMP REQUIRED)
find_package(CURL)
find_package(ZLIB)

# --- 复用你的实现：做一个 support 静态库 ---
add_library(halign4_test_support STATIC
        "${HALIGN4_ROOT}/src/utils/preprocess.cpp"
        "${HALIGN4_ROOT}/src/utils/seq_io.cpp"
        "${HALIGN4_ROOT}/src/utils/cmd.cpp"
        "${HALIGN4_ROOT}/src/utils/file_io.cpp"
        "${HALIGN4_ROOT}/src/mash/xxhash.c"
        "${HALIGN4_ROOT}/src/mash/hash.cpp"
        "${HALIGN4_ROOT}/src/consensus/selector.cpp"
        "${HALIGN4_ROOT}/src/consensus/consensus.cpp"
        "${HALIGN4_ROOT}/src/seed/minimizer.cpp"
        "${HALIGN4_ROOT}/src/seed/anchor.cpp"
        "${HALIGN4_ROOT}/src/mash/mash.cpp"
        "${HALIGN4_ROOT}/src/align/ref_aligner.cpp"
        "${HALIGN4_ROOT}/src/align/ksw2_extz2_sse.c"
        "${HALIGN4_ROOT}/src/align/psw_gg3_sse.c"
        "${HALIGN4_ROOT}/src/align/cigar.cpp"
        "${HALIGN4_ROOT}/src/align/align.cpp"
)

target_include_directories(halign4_test_support PUBLIC
        "${HALIGN4_ROOT}/include"
        "${HALIGN4_ROOT}/extern/CLI11/include"
        "${HALIGN4_ROOT}/extern/spdlog/include"
        "${HALIGN4_ROOT}/extern/cereal/include"
        "${HALIGN4_ROOT}/extern/doctest/include"
        "${HALIGN4_ROOT}/extern/WFA2-lib"
)

# 说明：这里显式链接 WFA2 静态库，避免只暴露 include 目录却没有对应库，
#       否则使用 WFA2 接口的测试在链接阶段可能报未定义符号。
target_link_libraries(halign4_test_support PUBLIC
        OpenMP::OpenMP_CXX
        wfa2::wfa2cpp_static
)

# 在 Release 时对测试支持库也启用与主项目一致的激进优化
if (MSVC)
    target_compile_options(halign4_test_support PRIVATE
        $<$<CONFIG:Release>:/O2>
        $<$<CONFIG:Release>:/arch:AVX2>
        $<$<CONFIG:Release>:/GL>
    )
    target_link_options(halign4_test_support PRIVATE
        $<$<CONFIG:Release>:/LTCG>
    )
else()
    target_compile_options(halign4_test_support PRIVATE
        $<$<CONFIG:Release>:-O3>
        $<$<CONFIG:Release>:-march=native>
    )
endif()

if (CURL_FOUND)
    target_link_libraries(halign4_test_support PUBLIC CURL::libcurl)
    target_compile_definitions(halign4_test_support PUBLIC HALIGN4_HAVE_LIBCURL=1)
endif()

if (ZLIB_FOUND)
    target_link_libraries(halign4_test_support PUBLIC ZLIB::ZLIB)
    target_compile_definitions(halign4_test_support PUBLIC HALIGN4_HAVE_ZLIB=1)
endif()

target_compile_definitions(halign4_test_support PUBLIC
        $<$<CONFIG:Debug>:_DEBUG>
)

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

# 注意：本 CMakeLists 位于 test/ 目录下，这里需要返回到工程根目录再进入 extern/WFA2-lib。
# CMake 规定：add_subdirectory() 引入“非当前目录的子目录”时，必须显式指定 binary dir，
# 否则会报："not a subdirectory ... must be explicitly specified"。
#
# 这里把 WFA2-lib 的构建目录固定放到 test 的构建树里，避免污染主工程 build 目录，
# 也避免多处 add_subdirectory 相互冲突。
add_subdirectory(${HALIGN4_ROOT}/extern/WFA2-lib ${CMAKE_CURRENT_BINARY_DIR}/wfa2lib-build)

# 注意：这里不能直接链接到主工程 target `halign4`（在此 CMakeLists 中并未定义），
# 否则在配置阶段会报错。我们将 OpenMP + WFA2 静态库链接到测试支撑库，
# 让 halign4_tests 通过 halign4_test_support 间接获得这些依赖，保持测试逻辑不变。
target_link_libraries(halign4_test_support PRIVATE OpenMP::OpenMP_CXX wfa2::wfa2cpp_static)


add_executable(halign4_tests
        "${CMAKE_CURRENT_LIST_DIR}/doctest_main.cpp"
        "${CMAKE_CURRENT_LIST_DIR}/test_consensus.cpp"
        "${CMAKE_CURRENT_LIST_DIR}/test_read_fasta.cpp"
        "${CMAKE_CURRENT_LIST_DIR}/test_minimizer.cpp"
        "${CMAKE_CURRENT_LIST_DIR}/test_mash.cpp"
        "${CMAKE_CURRENT_LIST_DIR}/test_jaccard.cpp"
        "${CMAKE_CURRENT_LIST_DIR}/test_fasta_writer.cpp"
        "${CMAKE_CURRENT_LIST_DIR}/test_align.cpp"
        "${CMAKE_CURRENT_LIST_DIR}/test_anchor.cpp"
)

target_link_libraries(halign4_tests PRIVATE halign4_test_support)

# 对测试可执行文件也应用 Release 优化选项
if (MSVC)
    target_compile_options(halign4_tests PRIVATE
        $<$<CONFIG:Release>:/O2>
        $<$<CONFIG:Release>:/arch:AVX2>
        $<$<CONFIG:Release>:/GL>
    )
    target_link_options(halign4_tests PRIVATE
        $<$<CONFIG:Release>:/LTCG>
    )
else()
    target_compile_options(halign4_tests PRIVATE
        $<$<CONFIG:Release>:-O3>
        $<$<CONFIG:Release>:-march=native>
        $<$<CONFIG:Release>:-flto>
        $<$<CONFIG:Release>:-pipe>
    )
    target_link_options(halign4_tests PRIVATE
        $<$<CONFIG:Release>:-flto>
    )
endif()

# 建议显式启用（即使 include(CTest) 已经做了也不冲突）
enable_testing()

# 注册一个测试：运行这个可执行文件
# add_test(NAME halign4_tests COMMAND halign4_tests)

add_test(NAME consensus
        COMMAND halign4_tests --test-suite=consensus)

add_test(NAME read_fasta
        COMMAND halign4_tests --test-suite=read_fasta)


add_test(NAME minimizer
        COMMAND halign4_tests --test-suite=minimizer)

add_test(NAME mash
        COMMAND halign4_tests --test-suite=mash)

add_test(NAME jaccard
        COMMAND halign4_tests --test-suite=jaccard)

add_test(NAME write_fasta
        COMMAND halign4_tests --test-suite=write_fasta)

# Align 正确性测试（快速）
add_test(NAME align
        COMMAND halign4_tests --test-suite=align)

# Align 性能测试（耗时较长，包含高相似度序列测试）
add_test(NAME align_perf
        COMMAND halign4_tests --test-suite=align_perf)

# Anchor 正确性测试（快速）
add_test(NAME anchor
        COMMAND halign4_tests --test-suite=anchor)


