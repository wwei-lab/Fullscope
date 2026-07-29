cmake_minimum_required(VERSION 3.16)

project(
  Fullscope
  VERSION 1.2.1
  DESCRIPTION "Full-length spatial transcriptomics processing toolkit"
  LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH ON)

find_package(Threads REQUIRED)

# Prefer pkg-config because it carries HTSlib's transitive link flags.  Some
# Conda environments, however, contain a usable libhts.so while htslib.pc
# cannot be resolved (most commonly because zlib.pc is missing).  Fall back to
# locating the shared library directly so those environments remain usable.
find_package(PkgConfig QUIET)
if(PkgConfig_FOUND)
  pkg_check_modules(HTSLIB QUIET IMPORTED_TARGET htslib)
endif()

if(TARGET PkgConfig::HTSLIB)
  set(FULLSCOPE_HTSLIB_TARGET PkgConfig::HTSLIB)
else()
  find_path(
    HTSLIB_INCLUDE_DIR
    NAMES htslib/hts.h
    HINTS "$ENV{CONDA_PREFIX}/include"
  )
  find_library(
    HTSLIB_LIBRARY
    NAMES hts
    HINTS "$ENV{CONDA_PREFIX}/lib"
  )

  if(NOT HTSLIB_INCLUDE_DIR OR NOT HTSLIB_LIBRARY)
    message(FATAL_ERROR
      "HTSlib was not found. Install the core dependencies with:\n"
      "  mamba install -c conda-forge -c bioconda htslib zlib\n"
      "Then rerun install.sh in the activated Conda environment."
    )
  endif()

  add_library(fullscope_htslib UNKNOWN IMPORTED)
  set_target_properties(
    fullscope_htslib
    PROPERTIES
      IMPORTED_LOCATION "${HTSLIB_LIBRARY}"
      INTERFACE_INCLUDE_DIRECTORIES "${HTSLIB_INCLUDE_DIR}"
  )
  set(FULLSCOPE_HTSLIB_TARGET fullscope_htslib)
  message(WARNING
    "pkg-config could not resolve HTSlib; using ${HTSLIB_LIBRARY} directly. "
    "Installing zlib in this Conda environment is recommended."
  )
endif()

find_package(seqan3 CONFIG REQUIRED)
find_package(cereal CONFIG REQUIRED)

set(FULLSCOPE_SOURCES
  scripts/src/CID_extract.cpp
  scripts/src/CID_extract_nogene.cpp
  scripts/src/CID_index.cpp
  scripts/src/CID_index_precise.cpp
  scripts/src/CID_mapping.cpp
  scripts/src/CID_mapping_precise.cpp
  scripts/src/MurmurHash3.cpp
  scripts/src/edlib.cpp
  scripts/src/fastq_segment.cpp
  scripts/src/main.cpp
  scripts/src/process_stereo_seq_bam.cpp
  scripts/src/utils.cpp
)

add_executable(fullscope-core ${FULLSCOPE_SOURCES})
target_include_directories(fullscope-core PRIVATE scripts/include)
target_link_libraries(
  fullscope-core
  PRIVATE
    ${FULLSCOPE_HTSLIB_TARGET}
    Threads::Threads
    seqan3::seqan3
    cereal::cereal
)

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  target_compile_options(fullscope-core PRIVATE -O3 -Wall -Wextra -Wpedantic)
endif()

install(TARGETS fullscope-core RUNTIME DESTINATION libexec/fullscope)

install(
  PROGRAMS
    fullscope_toolkit/bin/fullscope-batch
    fullscope_toolkit/bin/fullscope
    fullscope_toolkit/bin/fullscope-make-testdata
    fullscope_toolkit/bin/fullscope-merge
    fullscope_toolkit/bin/fullscope-submit
  DESTINATION bin
)

# Remove user-facing wrappers retired by the unified `fullscope` command.
install(CODE [[
  file(REMOVE
    "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/fullscope-ont"
    "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/fullscope-segment"
  )
]])

install(
  PROGRAMS
    fullscope_toolkit/scripts/fullscope_batch_submit.sh
    fullscope_toolkit/scripts/fullscope_make_testdata.sh
    fullscope_toolkit/scripts/fullscope_ont_pipeline.sh
    fullscope_toolkit/scripts/fullscope_segment.sh
    fullscope_toolkit/scripts/run_fullscope_ont_pipeline.slurm
  DESTINATION share/fullscope/scripts
)

install(
  FILES
    fullscope_toolkit/scripts/bambu_process.R
    fullscope_toolkit/scripts/fs_merge_spatial_transcripts.R
  DESTINATION share/fullscope/scripts
)

install(
  FILES
    fullscope_toolkit/refdata/adapters.fa
    fullscope_toolkit/refdata/anchor.fa
  DESTINATION share/fullscope/refdata
)

install(
  FILES
    fullscope_toolkit/config/config.example.env
    fullscope_toolkit/config/samples.example.tsv
  DESTINATION share/fullscope/config
)

install(FILES fullscope_toolkit/VERSION DESTINATION share/fullscope)
