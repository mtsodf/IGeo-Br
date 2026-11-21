# Build do iGeo com CMake (gfortran) e testes de verificação com CTest.
#
#   cmake -S . -B build                         # Release (padrão)
#   cmake --build build -j8                     # build/iGeo
#   ctest --test-dir build --output-on-failure  # testes de verificação
#
#   cmake -S . -B build-debug -DCMAKE_BUILD_TYPE=Debug   # -O0, -fcheck=all
#   cmake --build build -t BuildGraphs          # pós-processador de gráficos
#
# Opções:
#   -DIGEO_NATIVE=OFF      não otimiza para a CPU desta máquina
#   -DIGEO_TESTS=OFF       não registra os testes
#   -DBLA_VENDOR=OpenBLAS  escolhe a implementação de BLAS/LAPACK
#   -DIGEO_PYTHON=...      Python dos testes (padrão: .venv da raiz, se existir,
#                          senão python3 do PATH; ver requirements.txt)
#
# Equivale ao iGeo_vrs2.0/Makefile. No Windows com o compilador Intel,
# continue usando Build.bat.

cmake_minimum_required(VERSION 3.21)
project(iGeo VERSION 2.0 LANGUAGES Fortran)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Debug ou Release" FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Debug Release)
endif()

if(NOT CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  message(FATAL_ERROR
    "Este CMakeLists.txt é para o gfortran (compilador atual: "
    "${CMAKE_Fortran_COMPILER_ID}). Os módulos de compat/ substituem o IFPORT e "
    "o MKL da Intel; com ifort/ifx use Build.bat.")
endif()

option(IGEO_NATIVE "Otimiza para a CPU desta máquina (-mcpu/-march=native)" ON)
option(IGEO_TESTS "Registra os testes de verificação no CTest" ON)

find_package(OpenMP REQUIRED COMPONENTS Fortran)
# BLAS/LAPACK (dswap, dgeev): no macOS encontra o framework Accelerate.
find_package(LAPACK REQUIRED)

set(SRC ${CMAKE_CURRENT_SOURCE_DIR}/iGeo_vrs2.0)

# ---------------------------------------------------------------------------
# Flags (mesmas do Makefile; lá estão explicadas uma a uma)
# ---------------------------------------------------------------------------
# -ffree-line-length-none     linhas com mais de 132 colunas
# -fallow-argument-mismatch   chamadas ao PARDISO com escalares e vetores
# -fdollar-ok                 nome SRT$INTEGER4 (IFPORT)
# -frecursive                 variáveis locais na pilha (OpenMP)
# -finit-local-zero/-derived  o código conta com variáveis locais zeradas
add_library(igeo_flags INTERFACE)
target_compile_options(igeo_flags INTERFACE
  -ffree-line-length-none -fallow-argument-mismatch -fdollar-ok
  -frecursive -finit-local-zero -finit-derived
  # Sem -ffpe-trap: o código depende de aritmética IEEE (1/0 = infinito).
  $<$<CONFIG:Debug>:-O0 -g -fcheck=all -fbacktrace>
  $<$<CONFIG:Release>:-O3 -g>)
target_link_libraries(igeo_flags INTERFACE OpenMP::OpenMP_Fortran LAPACK::LAPACK)

if(IGEO_NATIVE)
  include(CheckFortranCompilerFlag)
  check_fortran_compiler_flag(-mcpu=native IGEO_HAS_MCPU_NATIVE)    # arm64
  if(IGEO_HAS_MCPU_NATIVE AND CMAKE_SYSTEM_PROCESSOR MATCHES "arm64|aarch64")
    target_compile_options(igeo_flags INTERFACE $<$<CONFIG:Release>:-mcpu=native>)
  else()
    check_fortran_compiler_flag(-march=native IGEO_HAS_MARCH_NATIVE) # x86_64
    if(IGEO_HAS_MARCH_NATIVE)
      target_compile_options(igeo_flags INTERFACE $<$<CONFIG:Release>:-march=native>)
    endif()
  endif()
endif()

# ---------------------------------------------------------------------------
# Simulador
# ---------------------------------------------------------------------------
# Programas independentes de pós-processamento (não fazem parte do iGeo).
set(GRAPH_SRCS ${SRC}/InjeSim_BuildGaphs.f90 ${SRC}/InjeSim_BuildGaphs2.f90)

file(GLOB IGEO_SRCS CONFIGURE_DEPENDS ${SRC}/*.f90)
list(REMOVE_ITEM IGEO_SRCS ${GRAPH_SRCS})
# Substitutos do IFPORT e do Intel MKL.
list(APPEND IGEO_SRCS ${SRC}/compat/ifport.f90 ${SRC}/compat/mkl_compat.f90)

# O CMake ordena a compilação pelas dependências entre módulos (USE).
add_executable(iGeo ${IGEO_SRCS})
target_link_libraries(iGeo PRIVATE igeo_flags)
set_target_properties(iGeo PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod/iGeo
  RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

# Pós-processador de gráficos (BuildGraphs.bat usa InjeSim_BuildGaphs2.f90).
add_executable(BuildGraphs EXCLUDE_FROM_ALL ${SRC}/InjeSim_BuildGaphs2.f90)
target_link_libraries(BuildGraphs PRIVATE igeo_flags)
set_target_properties(BuildGraphs PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod/BuildGraphs
  RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

# ---------------------------------------------------------------------------
# Testes
# ---------------------------------------------------------------------------
if(IGEO_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif()
