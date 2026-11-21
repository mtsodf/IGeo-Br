+++++++++++++++++++++++++++++++++++++++++++++
iGeo - caso de verificação
Injeção em fratura de comprimento fixo (deformação plana)
+++++++++++++++++++++++++++++++++++++++++++++

Solução de referência
  Abertura (Sneddon 1946; Sneddon & Lowengrub 1969), fratura de meio
  comprimento L sob pressão líquida uniforme Dp, em deformação plana:
      w(x) = 4 (1-nu^2) Dp / E * sqrt(L^2 - x^2)      (abertura total)
  O modelo é 1/4 do domínio (x >= 0, y >= 0, simetria em x = 0 e y = 0),
  então o iGeo calcula uy = w/2 na face y = 0.
  A injeção é feita a pressão controlada no poço: com a fratura muito
  condutiva, a pressão dentro dela é uniforme e igual à do poço, e
  Dp = p_poço - Shmin. Volume de fluido na fratura (quarto de modelo):
      V = H * Integral_0^L uy dx = pi (1-nu^2) Dp L^2 H / (2 E)
  (H = altura da fratura = espessura da camada)

  Obs. 1: a injeção a vazão controlada não é usada porque o acoplamento
  sequencial fluxo -> mecânica do iGeo só considera a complacência da
  fratura no passo seguinte; em rocha impermeável a pressão diverge.
  Obs. 2: os testes também rodam uma variante com a pressão mantida por
  vários passos, para verificar que a abertura não deriva.

Hipóteses representadas
  - deformação plana: uma camada (NZ = 1) com uz = 0 no topo e na base;
  - rocha elástica: Kg = Kb = E/(3(1-2nu)) => coeficiente de Biot ~ 0;
  - rocha praticamente impermeável (1e-5 mD): perda de fluido desprezível;
  - pressão inicial igual a Shmin: pressão líquida inicial nula;
  - fratura fixa: sem #CDM_FRACTURE, não há propagação;
  - abertura inicial de 1 mm só para a condutividade (não entra no volume),
    de modo que a pressão dentro da fratura seja praticamente uniforme.

Parâmetros: E = 1.0E7 kPa, nu = 0.25, L = 5 m, H = 1 m, Shmin = 10000 kPa,
pressão no poço 11000 kPa (Dp = 1000 kPa). Malha: dx = 0.25 m (20 elementos
na fratura), domínio de 20L x 20L. Resultado esperado (Sneddon):
uy(0) = 2 (1-nu^2) Dp L / E = 0.9375 mm e V = 3.6816E-3 m3.

*************************************************************************
* SOLVERS
*************************************************************************
#NTHREADS
4

#FLOWFLAG
#STRUCTURALFLAG

#STRUCT_ViscoplasticSolver
100000              <- max. iterações do solver linear
1.0E-9              <- tolerância do solver linear
25                  <- max. passos viscoplásticos
1.0E-2              <- tolerância do solver viscoplástico

*************************************************************************
* FLUIDO
*************************************************************************
#FLUIDINITIALIZATION
-500.0   10000.0    <- OWdepth(m)  pressão(kPa)
0.0      0.0        <- oilGrad(kPa/m)  waterGrad(kPa/m)

#WATERPHASEPROPS
1.0                 <- water_bwI
10000.0             <- water_pref (kPa)
1.0E-6              <- density_pref (desprezível: pressão uniforme na vertical)
1.0E-7              <- water_cw (1/kPa)
1.0                 <- water_viscI (cP)

*************************************************************************
* FRATURA
*************************************************************************
#PREFRACTURE
1                   <- número de completações
1  1  5.0  1.0E-3   <- K1, K2, Xlength (m), abertura inicial (m)

*************************************************************************
* WORKFLOW: INJEÇÃO
*************************************************************************
#WF_WaterInj

#WELLPRESSURECONSTRAINT
-999.0  -1001.0     <- z1, z2
>                   <- tempo (min), pressão no poço (kPa), Zref (m)
1.0    11000.0  -1000.0
<

*************************************************************************
* MALHA E PROPRIEDADES
*************************************************************************
Colunas: Nome Topo(m) POR Kh(mD) Kv(mD) Eh(kPa) Ev(kPa) niH niV Kg(GPa) Sv Shmax Shmin (kPa)
#NWB_STRUCTUREDMESH
>
SneddonInj          <- nome do poço (prefixo das saídas)
BOX                 <- esquema de malha
400  45  1  1.15    <- NX, NY, NZ, razão de progressão em Y
100.0  100.0        <- LX, LY (m)
#TABLEofProperties
1  -1000.0  0.10  1.0E-5  1.0E-5  1.0E7  1.0E7  0.25  0.25  6.666667  15000.0  12000.0  10000.0
<
