# iGeo

**iGeo** é um simulador de geomecânica computacional por elementos finitos, escrito em Fortran 90, para análise de poços de petróleo e de seus reservatórios. Ele acopla, em malhas 3D de hexaedros:

- **mecânica dos sólidos**: elasticidade, elastoplasticidade, viscoplasticidade e fluência;
- **escoamento em meio poroso**: fluxo transiente e permanente;
- **transporte de calor**: acoplamento termo-hidro-mecânico (THM);
- **propagação de fraturas hidráulicas**: mecânica do dano contínuo.

O código foi desenvolvido por **Luis Carlos de Sousa Junior** no Programa de Pós-Graduação em Geotecnia da **Universidade de Brasília (UnB)**, com orientação de **Manoel Porfírio Cordão Neto**, e serviu de base para a tese de doutorado e para artigos sobre injeção de água e fraturamento hidráulico.

---

## Principais funcionalidades

| Área | Recursos |
|---|---|
| **Malhas** | Malhas estruturadas próximas ao poço (*near wellbore*) em esquema `BOX` (cartesiano) ou `RAD` (radial), poços verticais e horizontais, elementos hexaédricos de 8 nós |
| **Modelos constitutivos** | Elástico (isotrópico, transversalmente isotrópico, ortotrópico), Mohr-Coulomb, Cam-Clay Modificado, modelo hiperbólico para areias, dano por tração e por cisalhamento, fluência de sal (*salt creep*) |
| **Solvers estruturais** | Viscoplástico implícito (padrão, `Viscoplastic`) e dinâmico explícito multifísico (`DynamicFEM` / `#MEDSolver`) |
| **Escoamento** | Fluxo monofásico ligeiramente compressível, condições de contorno de pressão ou de vazão no poço, injetividade, evolução de permeabilidade com a deformação plástica e dano à formação |
| **Térmico** | Solução numérica do campo de temperatura ou zona térmica analítica de Perkins & Gonzales |
| **Fratura hidráulica** | Fraturas iniciais (*pre-fracture*), propagação por dano à tração (CDM, energia de fratura `Gc`), transmissibilidade na fratura, fechamento de fratura |
| **Perfuração** | Escavação do poço, peso de fluido de perfuração (*mud weight*), estabilidade de poço |
| **Saída** | Arquivos VTK por passo de tempo (para ParaView/VisIt), relatórios texto e arquivo binário `.femOutput` com todo o histórico |
| **Paralelismo** | OpenMP e Intel MKL (PARDISO, `mkl_dcsrgemv`), além de gradiente conjugado pré-condicionado (PCG) com matrizes CSR |

---

## Estrutura do repositório

```
IGeo-Br/
├── iGeo_vrs2.0/          # Código-fonte (Fortran 90) e scripts de build/execução (Windows .bat)
├── Simulation_cases/     # Arquivos de entrada (.geo) de casos de teste e validação
├── tests/                # Testes de verificação (Python + CTest)
├── tools/                # Pós-processamento (perfil da abertura da fratura)
├── docs/                 # Relatórios (verificação de Sneddon)
├── CMakeLists.txt        # Build com CMake (gfortran) e testes com CTest
├── requirements.txt      # Dependências Python (testes e pós-processamento)
└── LICENSE.md            # Licença MIT
```

### Código-fonte (`iGeo_vrs2.0/`)

**Programa principal e fluxos de trabalho (*workflows*)**

| Arquivo | Descrição |
|---|---|
| `IGeo.f90` | Programa principal: lê os argumentos, detecta o tipo de modelo e despacha para o fluxo adequado |
| `GeomechWell_main_these.f90` | Fluxo `WellModel0`: monta o modelo near-wellbore e executa `#WF_WaterInj` ou `#WF_MEDSolver` |
| `WellStability.f90` | Fluxo `WellModel1`: estabilidade de poço e perfuração (`#WELLBOREGRID`) |
| `WF_WaterInjection.f90` | Injeção de água acoplada com geomecânica e propagação de fratura (`#WF_WaterInj`) |
| `WF_MEDSolver.f90` | Problema poromecânico com o solver dinâmico explícito multifísico (`#WF_MEDSolver`) |
| `WF_Excavation.f90` | Escavação/perfuração do poço (`#WF_Excavation`) |

**Definição do problema e montagem do modelo**

| Arquivo | Descrição |
|---|---|
| `ProblemDefinition.f90` | Leitura do arquivo `.geo` (palavras-chave) e estrutura `Problem` |
| `NearWellboreModel.f90` | Geração da malha estruturada em torno do poço (BOX/RAD) |
| `WellModelBuilder.f90` | Montagem do modelo de poço para estabilidade/perfuração |
| `DiscretizedModel.f90` | Modelo geomecânico discreto (nós, elementos, DOFs, propriedades e estado) |
| `MeshAuxiliares.f90`, `MeshIO.f90` | Utilitários de malha |
| `MainConstants.f90` | Constantes globais (dimensão, nós e DOFs por elemento etc.) |

**Mecânica dos sólidos e cargas**

| Arquivo | Descrição |
|---|---|
| `ConstitutiveGeoModels.f90` | Modelos constitutivos (MC, MCC, hiperbólico, dano, elástico) |
| `ViscoplasticExplicitySolver.f90` | Solver viscoplástico (baseado nas rotinas de Smith & Griffiths, *Programming the Finite Element Method*) com PCG |
| `DynamicSolver.f90` | Solver dinâmico explícito (*MEDSolver*) |
| `CreepMod.f90` | Fluência de sal |
| `Damage.f90` | Atualização de dano |
| `StiffnessMatrixCSR.f90`, `ElementStiffness.f90` | Matriz de rigidez em formato CSR (isotrópica, TI e ortotrópica) |
| `InternalForceCalculator.f90`, `InterfaceMechanics.f90` | Forças internas e condições de apoio |
| `GravityLoad.f90`, `PorepressureLoad.f90`, `StressLoad.f90`, `ThermalLoad.f90`, `UniformLoad.f90` | Vetores de carga |
| `StructuralModel.f90` | Leitura dos parâmetros do solver estrutural |

**Escoamento, calor e fratura**

| Arquivo | Descrição |
|---|---|
| `FluidFlowSolver.f90` | Escoamento em meio poroso e no tubo do poço, condições de contorno, transmissibilidades |
| `PCG_Flow.f90` | Solver PCG específico para o escoamento |
| `Thermal.f90` | Modelo térmico (numérico e Perkins & Gonzales) |
| `HydraulicFractureModel.f90` | Modelo de fratura hidráulica e propagação por dano à tração |

**Bibliotecas numéricas e E/S**

| Arquivo | Descrição |
|---|---|
| `FEMlib.f90`, `NumericalUtilities.f90`, `NumericalLib.f90`, `SolidMechanicsLib.f90` | Funções de forma, pontos de Gauss, matriz B, invariantes e autovalores de tensores, solver PARDISO |
| `ParallelRoutines.f90` | Operações vetoriais paralelas (OpenMP) |
| `IO.f90` | Mensagens, busca de palavras-chave e leitura de tabelas |
| `Report.f90`, `VTK.f90`, `VTKgeneral.f90` | Relatórios e exportação VTK |

**Pós-processamento**

| Arquivo | Descrição |
|---|---|
| `InjeSim_BuildGaphs.f90`, `InjeSim_BuildGaphs2.f90` | Programas independentes que leem o `.femOutput` e extraem perfis (p, p', q, deslocamento) ao longo de uma linha da malha para `OUTPUT.txt` |

---

## Requisitos

Há dois caminhos de build:

- **Windows, compilador Intel** (`Build.bat`): Intel Fortran (`ifort` ou `ifx`), Intel MKL e OpenMP.
- **macOS / Linux, gfortran** (`Makefile` ou `CMakeLists.txt`): GNU Fortran com OpenMP e BLAS/LAPACK; CMake ≥ 3.21 para o build com CMake. No macOS, o BLAS/LAPACK vem do framework Accelerate, que já faz parte do sistema. Com Homebrew: `brew install gcc`.

Os testes de verificação e as ferramentas de pós-processamento usam Python 3 com numpy (≥ 2.0) e matplotlib, listados em `requirements.txt`. Crie o ambiente na raiz do repositório:

```bash
python3 -m venv .venv
.venv/bin/pip install -r requirements.txt
```

O CMake e os Makefiles usam `.venv/bin/python3` quando esse ambiente existe; sem ele, usam o `python3` do PATH.

O Intel ifx não existe para macOS, então em Mac o único caminho é o gfortran.

Para visualizar os resultados, use um leitor de VTK legado, como o [ParaView](https://www.paraview.org/).

---

## Compilação

No *Intel oneAPI Command Prompt* (Windows), dentro de `iGeo_vrs2.0/`:

```bat
Build.bat     :: build otimizado: /O3 /Qm64 /Qpar /Qopenmp /Qmkl:parallel
Build0.bat    :: build sem otimização; antes apaga *.exe, *.mod, *.obj e *.log
```

Os dois scripts geram `iGeo.exe`. **A ordem dos arquivos na linha de compilação importa**, porque os módulos Fortran precisam ser compilados antes de quem os usa. Use a ordem de `Build.bat` como referência.

O utilitário de gráficos é compilado e executado com `BuildGraphs.bat`, que gera `BuildGaphs.exe`.

### gfortran (macOS / Linux) com `make`

Há um `Makefile` em `iGeo_vrs2.0/` para o **gfortran**. Ele foi testado com o GNU Fortran 16 (Homebrew) em macOS arm64.

```bash
cd iGeo_vrs2.0
make -j8                 # build/iGeo (otimizado: -O3 -mcpu=native, OpenMP)
make -j8 MODE=debug      # build/iGeo-debug (-O0, -fcheck=all, backtrace)
make graphs              # build/BuildGraphs (pós-processador)
make run CASE=LinearFlowTest.geo   # roda um caso de ../Simulation_cases
make clean
```

As dependências entre módulos são geradas automaticamente a partir das instruções `USE` (`mkdeps.awk`), então a ordem dos arquivos não importa e a compilação paralela funciona. Mudanças no próprio `Makefile` recompilam tudo.

O código foi escrito para o compilador Intel. Para compilar com gfortran sem alterar o build do Windows, `compat/` traz substitutos das dependências da Intel:

| Dependência Intel | Substituto no build gfortran |
|---|---|
| Módulo `IFPORT` (`SLEEP`, `IARGC`, `GETARG`, `SORTQQ`) | `compat/ifport.f90`: as três primeiras já existem no gfortran; `SORTQQ` foi implementada (heapsort) |
| MKL `mkl_dcsrgemv`, `ddoti`, `vdmul` | `compat/mkl_compat.f90`: implementações equivalentes (com OpenMP) |
| MKL PARDISO | sem substituto; a única rotina que o chama (`LinearSolverPAR`) não é usada |
| BLAS/LAPACK (`dswap`, `dgeev`) | framework Accelerate no macOS; `-llapack -lblas` em outros sistemas (`make LAPACK_LIBS=-lopenblas`, por exemplo) |

Os arquivos `compat/`, `Makefile` e `mkdeps.awk` não são usados pelo `Build.bat`.

### gfortran (macOS / Linux) com CMake

O `CMakeLists.txt` da raiz compila os mesmos fontes, com os mesmos flags do `Makefile`, e registra os testes de verificação no CTest. O CMake resolve sozinho a ordem dos módulos.

```bash
cmake -S . -B build                          # Release: -O3 -mcpu=native (ou -march=native), OpenMP
cmake --build build -j8                      # build/iGeo
cmake --build build -t BuildGraphs           # build/BuildGraphs (pós-processador)

cmake -S . -B build-debug -DCMAKE_BUILD_TYPE=Debug   # -O0, -fcheck=all, backtrace
cmake --build build-debug -j8
```

Opções: `-DIGEO_NATIVE=OFF` (não otimiza para a CPU da máquina), `-DBLA_VENDOR=OpenBLAS` (escolhe o BLAS/LAPACK; no macOS o padrão é o Accelerate), `-DIGEO_TESTS=OFF` (não registra os testes) e `-DIGEO_PYTHON=/caminho/python3` (interpretador dos testes; o padrão é o `.venv` da raiz, se existir, senão o `python3` do PATH).

---

## Execução

```bat
iGeo.exe <pasta_do_caso> <arquivo.geo>
```

Exemplo:

```bat
iGeo.exe "C:\casos\Simulation_cases" "MLS3B_injection.geo"
```

> ⚠️ Passe **sempre os dois argumentos**. Com menos de dois, o programa entra em modo de teste e usa um caminho e um arquivo fixos no código (`IGeo.f90`), que só existem na máquina do autor.

Os scripts `RunTests.bat`, `RunCase.bat` e `Run_USandInj.bat` mostram exemplos de execução, mas usam caminhos absolutos da máquina do autor e precisam ser ajustados antes do uso.

### Como o programa decide o que rodar

1. **Tipo de modelo**, definido pela palavra-chave de malha encontrada no `.geo`:
   - `#NWB_STRUCTUREDMESH` → `WellModel0`, o modelo near-wellbore (injeção e fratura);
   - `#WELLBOREGRID` → `WellModel1`, estabilidade de poço e perfuração.
2. **Workflow**, definido pela palavra-chave `#WF_...`:
   - `#WF_WaterInj`: injeção de água com fluxo, geomecânica, temperatura e fratura;
   - `#WF_MEDSolver`: solver dinâmico explícito multifísico;
   - `#WF_Excavation`: escavação do poço.
3. **Físicas ativas**, ligadas pelas flags `#FLOWFLAG`, `#STRUCTURALFLAG` e `#THERMALFLAG`.

### Arquivos de saída

Os arquivos são gravados na pasta do caso:

| Arquivo | Conteúdo |
|---|---|
| `<NomeDoPoço>.vtk`, `<NomeDoPoço>_<passo>.vtk` | Malha, dados estáticos e resultados por passo (pressão, deslocamento, tensões, deformações plásticas, abertura de fratura etc.) |
| `<NomeDoPoço>.REPORT.txt` | Relatório texto por passo |
| `<arquivo.geo>.femOutput` | Binário (Fortran `UNFORMATTED`) com malha e histórico completo, lido pelos utilitários `InjeSim_BuildGaphs*` |
| `<arquivo.geo>.log` | Log da execução, gravado no diretório corrente |

---

## Formato do arquivo de entrada (`.geo`)

O `.geo` é um arquivo texto livre, organizado por **palavras-chave** que começam com `#`:

- O leitor procura cada palavra-chave como **primeiro token da linha** e lê os dados das linhas seguintes. Texto fora das palavras-chave é ignorado e pode servir como comentário.
- Para **desativar** uma palavra-chave, coloque um `C` na frente dela, como em `C #THERMALFLAG`.
- Tabelas de comprimento variável ficam entre `>` e `<`.
- Unidades usuais: kPa, m, mD, °C, minutos e m³/d. Os comentários `<- ...` ao lado dos valores indicam o significado de cada coluna.

Exemplo resumido (baseado em `MLS3B_injection.geo`):

```text
#NTHREADS
3

#FLOWFLAG
#STRUCTURALFLAG
C #THERMALFLAG                              (desativado)

#FLUIDINITIALIZATION
-2300.70   23202.55                         <- OWdepth(m)  Pressure(kPa)
10.085     10.085                           <- oilGrad  waterGrad (kPa/m)

#CDM_FRACTURE
1
1  79  100.0  1.0E-1  455.0                 <- K1 K2 T0(kPa) Max.YYdilation Gc(J)

#WF_WaterInj
#WELLRATECONSTRAINT
-2497.64  -2509.99                          <- z1, z2 do intervalo injetor
>                                           <- tempo(min)  vazão(m³/d)  Zref(m)
0.1    1000  -2485.3
120    2000  -2485.3
<

#NWB_STRUCTUREDMESH
>
MLS3B_lower_A                               <- nome do poço
BOX                                         <- esquema de malha: BOX ou RAD
150  25  79  1.20                           <- Nrad, Ntheta, Nz, razão de progressão
150.0  150.0                                <- LX, LY
#TABLEofProperties
1  -2450.02  0.35  0  0  6760000  6760000  0.408  0.408  38  39884  39486  34187
...                                         <- camada, Topo, POR, Kh, Kv, Eh, Ev, niH, niV, Kg, Sv, SHmax, Shmin
<
```

Algumas palavras-chave importantes:

| Grupo | Palavras-chave |
|---|---|
| Geral | `#WORKPATH`, `#NTHREADS`, `#TIMESTEPS`, `#REPORT`, `#OUTPUT_MONITOR` |
| Físicas | `#FLOWFLAG`, `#STRUCTURALFLAG`, `#THERMALFLAG` |
| Solver estrutural | `#STRUCT_SOLVER`, `#STRUCT_LSOLVER`, `#STRUCT_ViscoplasticSolver`, `#STRUCT_DynaSOLVER`, `#MEDSolver`, `#CRITICALTIMESTEPFACTOR` |
| Malha e propriedades | `#NWB_STRUCTUREDMESH`, `#WELLBOREGRID`, `#TABLEofProperties`, `#TRASVERSEISOTROPIC`, `#EDITprops`, `#ELASTOPLASTIC_MODEL`, `#SALTSSCREEPMODEL` |
| Fluido | `#FLUIDINITIALIZATION`, `#WATERPHASEPROPS`, `#OILPHASEPROPS`, `#TRANSIENTFLOW`, `#OUTFLOW`, `#NoFlowBC`, `#PRESSUREBC` |
| Poço e injeção | `#WELLRATECONSTRAINT`, `#WELLPRESSURECONSTRAINT`, `#WATERINJECTIONRATE`, `#WATERINJECTIONPRESSURE`, `#WATERINJECTIONSTAGES`, `#VERTICALWELLPIPE`, `#HorizontallWell`, `#ProductionCompletion`, `#PERFORATIONDEPTH` |
| Térmico | `#InitialTemperature`, `#ThermalModel`, `#WaterInjectionTemperature`, `#PerkinsGonzales`, `#TemperatureVariation` |
| Fratura e dano | `#PREFRACTURE`, `#CDM_FRACTURE`, `#DamageShearStrain`, `#DamageLiquefaction`, `#PERMENHANCEMENT`, `#FORMATIONFLUIDDAMAGE`, `#FRAC_CLOUSURE` |
| Perfuração | `#MUDWEIGHTppg`, `#MudHistoryPPG`, `#DrillingInterval`, `#ProductionCasingDepth`, `#SECTIONREPORT` |

A referência definitiva de cada bloco é a rotina que o lê, em `ProblemDefinition.f90` e nos módulos `WF_*`.

---

## Testes de verificação

`tests/` contém casos com solução analítica conhecida, executados pelo iGeo e comparados automaticamente com a solução (Python 3 + numpy, `unittest`):

Com CMake, depois de compilar (os testes usam o executável do próprio build):

```bash
ctest --test-dir build --output-on-failure -j8         # todos os testes, ~20 s
ctest --test-dir build -L sneddon -V                   # só os de Sneddon, com a saída
ctest --test-dir build-debug --output-on-failure -j8   # com o build de debug, ~55 s
```

Com o `Makefile`:

```bash
make -C iGeo_vrs2.0 test               # compila e roda os testes
make -C iGeo_vrs2.0 test MODE=debug    # com o build de debug (verificação de limites)
```

No CTest, cada classe de teste é um teste: `sneddon.solucao_analitica` compara o caso com Sneddon, e `sneddon.convergencia_e_pressao_mantida` verifica a convergência com a malha e a abertura com a pressão mantida. Um teste pulado (executável não encontrado) conta como falha. Os casos usam 4 threads, e o CTest reserva 4 processadores por teste com `-j`.

| Caso | Solução de referência | O que é verificado |
|---|---|---|
| `tests/sneddon_injection/` | Fratura pressurizada em deformação plana (Sneddon, 1946; Sneddon & Lowengrub, 1969): `w(x) = 4(1−ν²)Δp/E·√(L²−x²)` | Injeção a pressão controlada numa fratura fixa (L = 5 m) em rocha elástica impermeável: pressão uniforme na fratura, abertura no centro (±2 %), perfil (erro L2 < 3 %), volume (±3 %), face fechada fora da fratura, redução do erro com o refino da malha e abertura constante com a pressão mantida |

Com a malha do teste (dx = 0,25 m, domínio 20L), o erro é de −1,4 % no centro, 1,9 % no perfil e −2,4 % no volume. Ele diminui com o refino da malha e com o tamanho do domínio, porque vem dos elementos lineares e do contorno fixo.

Esse caso levou a duas correções no carregamento da fratura, descritas em detalhe, com figuras, no relatório [`docs/relatorio_sneddon/`](docs/relatorio_sneddon/README.md) (versão em PDF: `make -C docs/relatorio_sneddon`):

- `ViscoplasticExplicitySolver.f90`: atribuição a um vetor ainda não alocado (falha de segmentação com gfortran).
- `HydraulicFractureModel.f90` (`UpdateHydraulicFracture`): antes, a cada passo o código reaplicava como carga na face o resíduo `p + σyy`, com σyy médio do elemento. No elemento da ponta a tensão é singular e o resíduo nunca se anulava, então a abertura crescia de passo em passo mesmo com a pressão constante (+10 % no 6º passo). Agora a primeira carga de cada elemento continua sendo `p + σyy`, que libera a tensão da face, e as seguintes aplicam só a variação da pressão do fluido.

Limitação que permanece: **injeção a vazão controlada em rocha impermeável não é estável.** O acoplamento sequencial fluxo → mecânica só considera a complacência da fratura no passo seguinte, e a pressão diverge. Por isso o caso usa pressão controlada.

---

## Pós-processamento: perfil da abertura da fratura

`tools/perfil_abertura_fratura.py` lê os VTK de passo (`<NomeDoPoço>_<passo>.vtk`), identifica o plano da fratura (a primeira fileira de elementos da malha BOX ou RAD) e gera uma figura com o perfil vertical da abertura junto ao poço, o perfil ao longo do comprimento e o mapa da abertura no plano da fratura. Também imprime um resumo por passo (abertura máxima, altura e comprimento abertos) e pode exportar os dados em CSV.

```bash
source .venv/bin/activate      # ambiente criado com requirements.txt (ver Requisitos)
python3 tools/perfil_abertura_fratura.py Simulation_cases --poco LinearFlow
python3 tools/perfil_abertura_fratura.py caso/ --z -2480 --passos 1-10 --csv abertura.csv
```

A figura é salva por padrão como `<NomeDoPoço>_perfil_abertura.png`, ao lado dos VTK. Use `--help` para ver todas as opções.

---

## Casos de simulação (`Simulation_cases/`)

| Arquivo | Descrição |
|---|---|
| `LinearFlowTest.geo`, `LinearFlow.geo` | Fluxo linear transiente em modelo 2D BOX (só escoamento) |
| `ThermalTest_Rad.geo` | Fluxo com solução térmica em malha radial |
| `Validation2layersModel.geo`, `Validation2layersModel_1.geo` | Validação de propagação de fratura em modelo de 2 camadas (solução de Simonsen et al.) |
| `FracPropag.geo` | Propagação de fratura em modelo 3D BOX elástico |
| `MLS3B_injection*.geo`, `M3_injection_test.geo` | Injeção em reservatório multicamadas com fratura, dano e solver viscoplástico |
| `MLS3B.geo` | Escavação/estabilidade de poço (`#WELLBOREGRID` + `#WF_Excavation`) |
| `OffshoreSand_HorWell.geo` | Poço horizontal em areia offshore com solver dinâmico |
| `iGeoV1_c1.geo` | Teste básico em modelo multicamadas |
| `FracPropag_test.geo` | *(arquivo vazio)* |

> Os casos trazem em `#WORKPATH` caminhos da máquina do autor. Na prática, a pasta usada é a que você passa como primeiro argumento do `iGeo.exe`.

---

## Limitações conhecidas

- No build com gfortran não há solver PARDISO. Isso não afeta os workflows atuais, que usam os solvers iterativos (PCG).
- O código depende de variáveis locais começarem zeradas (o Makefile usa `-finit-local-zero -finit-derived` para garantir isso) e de aritmética IEEE com infinitos (por isso o modo debug não usa `-ffpe-trap`).
- Há caminhos absolutos fixos em `IGeo.f90` (modo de teste) e nos scripts `.bat`.
- Todos os arquivos estão em UTF-8 (sem BOM), mas a maioria usa terminações de linha CRLF (Windows).
- Vários workflows antigos (uniaxial, triaxial, Terzaghi, carga lateral etc.) estão comentados em `GeomechWell_main_these.f90` e não estão disponíveis nesta versão.
- O repositório não tem testes automatizados. Os arquivos em `Simulation_cases/` servem como casos de verificação manual.

---

## Autores

- **Luis Carlos de Sousa Junior**, desenvolvedor (UnB, Pós-Graduação em Geotecnia)
- **Manoel Porfírio Cordão Neto**, orientador

## Licença

Distribuído sob a licença **MIT**. Veja [LICENSE.md](LICENSE.md).
