# Makefile do iGeo para gfortran (testado com GNU Fortran 16, macOS arm64).
#
# Alvos principais:
#   make              compila o simulador (build/iGeo)
#   make graphs       compila o pós-processador de gráficos (build/BuildGraphs)
#   make run CASE=arquivo.geo [CASEDIR=../Simulation_cases]
#   make test         compila e roda os testes de verificação (../tests)
#   make clean        remove o diretório build/
#
# Opções:
#   make MODE=debug   sem otimização, com verificações de execução
#                     (gera build/iGeo-debug)
#   make FC=gfortran-16 ...   escolhe outro executável do gfortran
#   make -j8          compilação paralela (as dependências entre módulos
#                     são geradas automaticamente a partir dos USE)
#
# No Windows com o compilador Intel, continue usando Build.bat.

FC      ?= gfortran
MODE    ?= release

BUILDDIR := build
OBJDIR   := $(BUILDDIR)/obj-$(MODE)
MODDIR   := $(BUILDDIR)/mod-$(MODE)
BINDIR   := $(BUILDDIR)

# O padrão do make define FC=f77; força o gfortran nesse caso.
ifeq ($(origin FC),default)
  FC := gfortran
endif

# ---------------------------------------------------------------------------
# Flags
# ---------------------------------------------------------------------------
# -ffree-line-length-none : o código tem linhas com mais de 132 colunas
#                           (o compilador Intel aceita sem reclamar).
# -fallow-argument-mismatch : chamadas a rotinas externas (PARDISO) com
#                           escalares e vetores no mesmo argumento.
# -fdollar-ok              : permite o nome SRT$INTEGER4 (IFPORT da Intel).
# -frecursive              : variáveis locais na pilha, necessário com OpenMP.
# -finit-local-zero -finit-derived : zera variáveis locais. O código tem
#                           variáveis locais usadas sem inicialização explícita
#                           e conta com elas valendo zero, como costuma
#                           acontecer no build Intel. (Os flags LOGICAL dos
#                           tipos derivados já têm valor inicial no código.)
COMMON_FLAGS := -fopenmp -ffree-line-length-none -fallow-argument-mismatch -fdollar-ok \
                -frecursive -finit-local-zero -finit-derived -J$(MODDIR) -I$(MODDIR)

ifeq ($(MODE),debug)
  # Sem -ffpe-trap: o código depende de aritmética IEEE (por exemplo,
  # 1/0 = infinito para o módulo volumétrico do fluido em elementos drenados).
  MODE_FLAGS := -O0 -g -fcheck=all -fbacktrace
else ifeq ($(MODE),release)
  MODE_FLAGS := -O3 -mcpu=native -g
else
  $(error MODE deve ser 'release' ou 'debug')
endif

FFLAGS  ?=
FCFLAGS := $(COMMON_FLAGS) $(MODE_FLAGS) $(FFLAGS)
LDFLAGS ?=

# BLAS/LAPACK (dswap, dgeev). No macOS usa o framework Accelerate, que já
# vem com o sistema. Em outros sistemas, por exemplo:
#   make LAPACK_LIBS="-lopenblas"
ifeq ($(shell uname -s),Darwin)
  LAPACK_LIBS ?= -framework Accelerate
else
  LAPACK_LIBS ?= -llapack -lblas
endif

LDLIBS  := $(LAPACK_LIBS) -fopenmp

# ---------------------------------------------------------------------------
# Fontes
# ---------------------------------------------------------------------------
# Programas independentes de pós-processamento (não fazem parte do iGeo).
GRAPH_SRCS := InjeSim_BuildGaphs.f90 InjeSim_BuildGaphs2.f90

# Substitutos do IFPORT e do Intel MKL (veja compat/).
COMPAT_SRCS := compat/ifport.f90 compat/mkl_compat.f90

SRCS := $(COMPAT_SRCS) $(filter-out $(GRAPH_SRCS),$(wildcard *.f90))
OBJS := $(addprefix $(OBJDIR)/,$(notdir $(SRCS:.f90=.o)))

# O modo debug gera executáveis com sufixo, para não sobrescrever o release.
SUFFIX := $(if $(filter debug,$(MODE)),-debug,)
IGEO   := $(BINDIR)/iGeo$(SUFFIX)
GRAPHS := $(BINDIR)/BuildGraphs$(SUFFIX)
DEPS   := $(OBJDIR)/deps.mk

vpath %.f90 . compat

# ---------------------------------------------------------------------------
# Regras
# ---------------------------------------------------------------------------
.PHONY: all graphs run test clean help
.DEFAULT_GOAL := all

all: $(IGEO)

$(IGEO): $(OBJS) | $(BINDIR)
	$(FC) $(FCFLAGS) $(LDFLAGS) -o $@ $(OBJS) $(LDLIBS)

graphs: $(GRAPHS)

# BuildGraphs.bat usa InjeSim_BuildGaphs2.f90.
$(GRAPHS): InjeSim_BuildGaphs2.f90 | $(BINDIR) $(MODDIR)
	$(FC) $(FCFLAGS) $(LDFLAGS) -o $@ $< $(LDLIBS)

# Mudanças no Makefile (flags) recompilam tudo.
$(OBJDIR)/%.o: %.f90 Makefile | $(OBJDIR) $(MODDIR)
	$(FC) $(FCFLAGS) -c -o $@ $<

# Dependências entre módulos, regeneradas quando algum fonte muda.
$(DEPS): $(SRCS) mkdeps.awk | $(OBJDIR)
	awk -v objdir=$(OBJDIR) -f mkdeps.awk pass=1 $(SRCS) pass=2 $(SRCS) > $@

ifeq ($(filter clean help,$(MAKECMDGOALS)),)
  -include $(DEPS)
endif

$(OBJDIR) $(MODDIR) $(BINDIR):
	mkdir -p $@

# Executa um caso: make run CASE=LinearFlowTest.geo
CASEDIR ?= ../Simulation_cases
run: $(IGEO)
	@test -n "$(CASE)" || { echo "Uso: make run CASE=arquivo.geo [CASEDIR=pasta]"; exit 1; }
	./$(IGEO) "$(CASEDIR)" "$(CASE)"

# Testes de verificação (Python 3 + numpy), usando o executável deste modo.
# Usa o ambiente ../.venv (requirements.txt) se existir; senão, python3 do PATH.
VENV_PYTHON := $(abspath ../.venv/bin/python3)
PYTHON ?= $(if $(wildcard $(VENV_PYTHON)),$(VENV_PYTHON),python3)
test: $(IGEO)
	IGEO_BIN=$(abspath $(IGEO)) $(PYTHON) -m unittest discover -s ../tests -v

clean:
	rm -rf $(BUILDDIR)

help:
	@sed -n '2,17p' Makefile | sed 's/^# \{0,1\}//'
