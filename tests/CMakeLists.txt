# Testes de verificação do iGeo (Python 3 + numpy), executados pelo CTest.
#
#   ctest --test-dir build --output-on-failure
#   ctest --test-dir build -L sneddon -V        # só os testes de Sneddon
#
# Cada teste roda o executável iGeo deste build (variável IGEO_BIN) numa
# pasta temporária e compara a saída com a solução de referência.

# Interpretador Python dos testes, em ordem de preferência:
#   1. -DIGEO_PYTHON=/caminho/python3, se indicado;
#   2. o ambiente .venv da raiz do repositório (ver requirements.txt);
#   3. o python3 do PATH, como o Makefile. (O FindPython3 não é usado porque
#      dá preferência a um ambiente conda definido, mesmo fora do PATH.)
set(IGEO_PYTHON "" CACHE FILEPATH
    "Interpretador Python 3 (com numpy) dos testes; vazio = .venv da raiz ou python3 do PATH")
if(IGEO_PYTHON)
  set(IGEO_PYTHON_EXE "${IGEO_PYTHON}")
else()
  find_program(IGEO_PYTHON_EXE NAMES python3 python
               PATHS "${PROJECT_SOURCE_DIR}/.venv/bin" "${PROJECT_SOURCE_DIR}/.venv/Scripts"
               NO_DEFAULT_PATH NO_CACHE)
  if(NOT IGEO_PYTHON_EXE)
    find_program(IGEO_PYTHON_EXE NAMES python3 python NO_CACHE)
  endif()
endif()
if(NOT IGEO_PYTHON_EXE)
  message(FATAL_ERROR "python3 não encontrado. Crie o ambiente com "
                      "'python3 -m venv .venv && .venv/bin/pip install -r requirements.txt', "
                      "use -DIGEO_PYTHON=... ou desligue os testes com -DIGEO_TESTS=OFF.")
endif()

execute_process(
  COMMAND "${IGEO_PYTHON_EXE}" -c "import numpy"
  RESULT_VARIABLE IGEO_SEM_NUMPY OUTPUT_QUIET ERROR_QUIET)
if(IGEO_SEM_NUMPY)
  message(WARNING "numpy não encontrado para ${IGEO_PYTHON_EXE}: os testes de "
                  "verificação vão falhar. Crie o ambiente da raiz com "
                  "'python3 -m venv .venv && .venv/bin/pip install -r requirements.txt'.")
else()
  message(STATUS "Testes de verificação com ${IGEO_PYTHON_EXE}")
endif()

# Registra uma classe de teste do unittest como um teste do CTest.
#   igeo_teste_python(<nome> <módulo> <classe> [TIMEOUT s] [LABELS ...])
function(igeo_teste_python nome modulo classe)
  cmake_parse_arguments(T "" "TIMEOUT" "LABELS" ${ARGN})
  if(NOT T_TIMEOUT)
    set(T_TIMEOUT 900)
  endif()
  add_test(NAME ${nome}
           COMMAND ${IGEO_PYTHON_EXE} -m unittest -v ${modulo}.${classe}
           WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
  set_tests_properties(${nome} PROPERTIES
    ENVIRONMENT "IGEO_BIN=$<TARGET_FILE:iGeo>;PYTHONDONTWRITEBYTECODE=1"
    TIMEOUT ${T_TIMEOUT}
    LABELS "verificacao;${T_LABELS}"
    # O .geo usa #NTHREADS 4: o CTest reserva 4 processadores com -j.
    PROCESSORS 4
    # O unittest sai com código 0 quando pula testes (executável não
    # encontrado); aqui isso é uma falha.
    FAIL_REGULAR_EXPRESSION "skipped=[0-9]")
endfunction()

# Fratura pressurizada de comprimento fixo: solução de Sneddon
# (tests/test_sneddon_injection.py, caso em tests/sneddon_injection/).
igeo_teste_python(sneddon.solucao_analitica test_sneddon_injection
                  TesteInjecaoSneddon LABELS sneddon)
igeo_teste_python(sneddon.convergencia_e_pressao_mantida test_sneddon_injection
                  TesteConvergenciaEPressaoMantida LABELS sneddon)
