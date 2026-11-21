+++++++++++++++++++++++++++++++++++++++++++++
iGeo - caso de verificação
Fratura PKN: altura fixa, propagação por injeção a vazão constante
+++++++++++++++++++++++++++++++++++++++++++++

Solução de referência (Nordgren, 1972; Perkins & Kern, 1961)
  Fratura vertical de altura H fixa, rocha elástica homogênea, fluido
  newtoniano, sem perda de fluido e sem tenacidade. Cada seção vertical
  está em deformação plana (seção elíptica). Com q = vazão de uma asa:
      L(t)    = 0.68 [G q^3 / ((1-nu) mu H^4)]^(1/5) t^(4/5)
      w0(t)   = 2.5  [(1-nu) mu q^2 / (G H)]^(1/5) t^(1/5)   (abertura no poço)
      pnet(t) = G w0 / ((1-nu) H)
  O modelo é 1/4 do domínio (simetria em x = 0 e y = 0): ele representa
  metade de uma asa (o lado y >= 0 da abertura), então a vazão no modelo
  é q/2 e o iGeo calcula uy = w/2 na face y = 0.

Hipóteses representadas
  - altura fixa: #CDM_FRACTURE só nas camadas 9 a 18 (H = 10 m, camadas
    de 1 m); a fratura não cresce para as camadas vizinhas;
  - rocha elástica e homogênea em todas as camadas (mesmos E, nu e
    tensões), com Kg = Kb => coeficiente de Biot ~ 0;
  - sem perda de fluido: permeabilidade de 1e-5 mD;
  - pressão inicial igual a Shmin: pressão líquida inicial nula;
  - fluido newtoniano de 100 cP, densidade desprezível (sem gravidade);
  - tenacidade desprezível: Gc = 1 J/m2 (K_IC ~ 0.1 MPa.m^0.5, contra uma
    pressão líquida de ~2.5 MPa);
  - pré-fratura de 2 m (um elemento) em toda a altura, com abertura
    inicial de 0.1 mm só para a condutividade.

Parâmetros: E = 1.0E7 kPa (G = 4 GPa), nu = 0.25, H = 10 m, mu = 100 cP,
q = 0.01 m3/s por asa (432 m3/d no modelo), Shmin = 10000 kPa.
Malha: dx = 2 m (domínio de 120 m), 20 elementos em y (60 m, progressão
1.25), 26 camadas em z (zona de 10 m e 8 camadas crescentes acima e abaixo).
Resultado esperado (Nordgren) em t = 120 s: L = 43.8 m, w0 = 4.66 mm,
pnet = 2.49 MPa, volume no modelo = 0.6 m3.

Obs.: o tratamento da ponta no iGeo (critério de propagação e zona seca à
frente do fluido) ainda não reproduz a ponta do PKN; veja
docs/relatorio_pkn/README.md. O teste verifica as partes que o iGeo deve
reproduzir (balanço de massa, seção elíptica, lei cúbica, contenção da
altura) e registra o desvio em relação a Nordgren como falha esperada.

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
1.0E-6              <- density_pref (desprezível: sem gravidade)
1.0E-7              <- water_cw (1/kPa)
100.0               <- water_viscI (cP)

*************************************************************************
* FRATURA
*************************************************************************
#CDM_FRACTURE
1                   <- número de regiões
9  18  0.0  1.0E-1  1.0    <- K1, K2, T0 (kPa), Max.YYdilation, Gc (J/m2)

#PREFRACTURE
1                   <- número de completações
9  18  2.0  1.0E-4  <- K1, K2, Xlength (m), abertura inicial (m)

*************************************************************************
* WORKFLOW: INJEÇÃO A VAZÃO CONSTANTE (passos de 10 s)
*************************************************************************
#WF_WaterInj

#WELLRATECONSTRAINT
-1036.000  -1045.000  <- z1, z2 (nós internos da zona: sem os da borda)
>                    <- tempo (min), vazão no modelo (m3/d), Zref (m)
0.166667  432.0  -1035.500
0.333333  432.0  -1035.500
0.500000  432.0  -1035.500
0.666667  432.0  -1035.500
0.833333  432.0  -1035.500
1.000000  432.0  -1035.500
1.166667  432.0  -1035.500
1.333333  432.0  -1035.500
1.500000  432.0  -1035.500
1.666667  432.0  -1035.500
1.833333  432.0  -1035.500
2.000000  432.0  -1035.500
<

*************************************************************************
* MALHA E PROPRIEDADES
*************************************************************************
Colunas: Nome Topo(m) POR Kh(mD) Kv(mD) Eh(kPa) Ev(kPa) niH niV Kg(GPa) Sv Shmax Shmin (kPa)
#NWB_STRUCTUREDMESH
>
PKN                 <- nome do poço (prefixo das saídas)
BOX                 <- esquema de malha
60  20  26  1.25    <- NX, NY, NZ, razão de progressão em Y
120.0  60.0         <- LX, LY (m)
#TABLEofProperties
1  -1000.000  0.10  1.0E-5  1.0E-5  1.0000E+07  1.0000E+07  0.25  0.25  6.666667  15000.0  12000.0  10000.0
2  -1010.000  0.10  1.0E-5  1.0E-5  1.0000E+07  1.0000E+07  0.25  0.25  6.666667  15000.0  12000.0  10000.0
3  -1018.000  0.10  1.0E-5  1.0E-5  1.0000E+07  1.0000E+07  0.25  0.25  6.666667  15000.0  12000.0  10000.0
4  -1024.000  0.10  1.0E-5  1.0E-5  1.0000E+07  1.0000E+07  0.25  0.25  6.666667  15000.0  12000.0  10000.0
5  -1028.000  0.10  1.0E-5  1.0E-5  1.0000E+07  1.0000E+07  0.25  0.25  6.666667  15000.0  12000.0  10000.0
6  -1031.000  0.10  1.0E-5  1.0E-5  1.0000E+07  1.0000E+07  0.25  0.25  6.666667  15000.0  12000.0  10000.0
7  -1033.000  0.10  1.0E-5  1.0E-5  1.0000E+07  1.0000E+07  0.25  0.25  6.666667  15000.0  12000.0  10000.0
8  -1034.500  0.10  1.0E-5  1.0E-5  1.0000E+07  1.0000E+07  0.25  0.25  6.666667  15000.0  12000.0  10000.0
9  -1035.500  0.10  1.0E-5  1.0E-5  1.0000E+07  1.0000E+07  0.25  0.25  6.666667  15000.0  12000.0  10000.0
10  -1036.500  0.10  1.0E-5  1.0E-5  1.0000E+07  1.0000E+07  0.25  0.25  6.666667  15000.0  12000.0  10000.0
11  -1037.500  0.10  1.0E-5  1.0E-5  1.0000E+07  1.0000E+07  0.25  0.25  6.666667  15000.0  12000.0  10000.0
12  -1038.500  0.10  1.0E-5  1.0E-5  1.0000E+07  1.0000E+07  0.25  0.25  6.666667  15000.0  12000.0  10000.0
13  -1039.500  0.10  1.0E-5  1.0E-5  1.0000E+07  1.0000E+07  0.25  0.25  6.666667  15000.0  12000.0  10000.0
14  -1040.500  0.10  1.0E-5  1.0E-5  1.0000E+07  1.0000E+07  0.25  0.25  6.666667  15000.0  12000.0  10000.0
15  -1041.500  0.10  1.0E-5  1.0E-5  1.0000E+07  1.0000E+07  0.25  0.25  6.666667  15000.0  12000.0  10000.0
16  -1042.500  0.10  1.0E-5  1.0E-5  1.0000E+07  1.0000E+07  0.25  0.25  6.666667  15000.0  12000.0  10000.0
17  -1043.500  0.10  1.0E-5  1.0E-5  1.0000E+07  1.0000E+07  0.25  0.25  6.666667  15000.0  12000.0  10000.0
18  -1044.500  0.10  1.0E-5  1.0E-5  1.0000E+07  1.0000E+07  0.25  0.25  6.666667  15000.0  12000.0  10000.0
19  -1045.500  0.10  1.0E-5  1.0E-5  1.0000E+07  1.0000E+07  0.25  0.25  6.666667  15000.0  12000.0  10000.0
20  -1046.500  0.10  1.0E-5  1.0E-5  1.0000E+07  1.0000E+07  0.25  0.25  6.666667  15000.0  12000.0  10000.0
21  -1048.000  0.10  1.0E-5  1.0E-5  1.0000E+07  1.0000E+07  0.25  0.25  6.666667  15000.0  12000.0  10000.0
22  -1050.000  0.10  1.0E-5  1.0E-5  1.0000E+07  1.0000E+07  0.25  0.25  6.666667  15000.0  12000.0  10000.0
23  -1053.000  0.10  1.0E-5  1.0E-5  1.0000E+07  1.0000E+07  0.25  0.25  6.666667  15000.0  12000.0  10000.0
24  -1057.000  0.10  1.0E-5  1.0E-5  1.0000E+07  1.0000E+07  0.25  0.25  6.666667  15000.0  12000.0  10000.0
25  -1063.000  0.10  1.0E-5  1.0E-5  1.0000E+07  1.0000E+07  0.25  0.25  6.666667  15000.0  12000.0  10000.0
26  -1071.000  0.10  1.0E-5  1.0E-5  1.0000E+07  1.0000E+07  0.25  0.25  6.666667  15000.0  12000.0  10000.0
<
