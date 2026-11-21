# Verificação do iGeo com a solução de Sneddon e correção do carregamento da fratura

## Resumo

- **Caso usado:** injeção, a pressão controlada, numa fratura de comprimento fixo, em deformação plana, numa rocha elástica impermeável. A solução analítica é a de Sneddon (1946) para uma fratura pressurizada.
- **Problema encontrado:** no primeiro passo de tempo o iGeo reproduzia a solução. A partir do segundo passo, a abertura crescia continuamente **mesmo com a pressão constante**: +11,7 % na abertura no centro e +22 % no volume depois de 6 passos.
- **Causa:** a cada passo, `UpdateHydraulicFracture` reaplicava como carga na face o resíduo `p + σyy`, com σyy igual à tensão **média** do elemento junto à face. No elemento da ponta a tensão é singular, essa média nunca chega a −p e o resíduo nunca zera.
- **Correção:** a primeira carga de cada elemento continua sendo `p + σyy`, que libera a tensão da face quando ela abre. As cargas seguintes aplicam só a **variação da pressão do fluido** desde a carga anterior.
- **Resultado:** a abertura fica constante enquanto a pressão é mantida. O erro em relação a Sneddon é de −1,4 % no centro e −2,4 % no volume, e diminui com o refino da malha e com o tamanho do domínio.
- **Efeito em fraturas que se propagam:** no caso `M3_injection_test`, a fratura passa a se propagar durante a injeção (43 eventos em 22 min, contra 1 antes), e a pressão cai em dente de serra a cada evento. Ainda falta validar a propagação contra uma solução de referência (por exemplo, o KGD).

---

## 1. O caso de referência

### 1.1 Solução de Sneddon

Uma fratura de meio comprimento *L*, num meio elástico infinito em deformação plana, sob pressão líquida uniforme Δp, tem abertura total elíptica (Sneddon, 1946; Sneddon & Lowengrub, 1969):

$$
w(x) = \frac{4\,(1-\nu^2)\,\Delta p}{E}\,\sqrt{L^2 - x^2}
$$

O iGeo modela 1/4 do domínio (x ≥ 0, y ≥ 0), com simetria em x = 0 e em y = 0. A grandeza que ele calcula é o deslocamento da face, igual à **meia abertura** `uy = w/2`. O volume de fluido na fratura, no quarto de modelo com altura *H*, é:

$$
V = H\int_0^L u_y\,dx = \frac{\pi\,(1-\nu^2)\,\Delta p\,L^2 H}{2E}
$$

### 1.2 Modelo no iGeo

Arquivo: [`tests/sneddon_injection/sneddon_injection.geo`](../../tests/sneddon_injection/sneddon_injection.geo).

![Esquema do caso](figuras/fig1_esquema.svg)

*Figura 1 — Quarto de modelo (desenho fora de escala). A fratura ocupa a face y = 0 entre o poço (x = 0) e x = L, e recebe a pressão do fluido. O resto da face y = 0 e os contornos têm deslocamento normal restrito.*

| Hipótese da solução | Como foi representada no iGeo |
|---|---|
| Deformação plana | Uma única camada (NZ = 1); o modelo BOX restringe `uz` no topo e na base |
| Rocha elástica | `Kg` igual ao módulo volumétrico drenado (`Kb = E/(3(1−2ν))`), então o coeficiente de Biot ≈ 0 |
| Sem perda de fluido para a rocha | Permeabilidade de 10⁻⁵ mD |
| Pressão líquida uniforme e conhecida | Pressão imposta no poço (`#WELLPRESSURECONSTRAINT`) e abertura inicial de 1 mm, que só entra na condutividade da fratura |
| Fratura fixa | Sem `#CDM_FRACTURE`, portanto sem propagação |
| Pressão líquida inicial nula | Pressão de poros inicial igual a Shmin (10 MPa) |
| Meio infinito | Domínio de 100 m × 100 m (20 L) |

Parâmetros: E = 10 GPa, ν = 0,25, L = 5 m, H = 1 m, pressão no poço de 11 MPa (Δp = 1 MPa). Resultado esperado: `uy(0) = 2(1−ν²)ΔpL/E = 0,9375 mm` e `V = 3,68×10⁻³ m³`.

![Malha junto à fratura](figuras/fig2_malha.svg)

*Figura 2 — Malha em planta junto à fratura: 400 × 45 elementos, dx = 0,25 m (20 elementos na fratura), com Y em progressão geométrica (o primeiro elemento tem cerca de 1 cm).*

**Por que pressão controlada, e não vazão.** O iGeo resolve o fluxo e a mecânica em sequência, e a complacência da fratura (o quanto ela abre por unidade de pressão) só entra no passo seguinte, como correção de volume. Em rocha impermeável, no primeiro passo de fluxo a fratura ainda tem volume zero e a pressão diverge (chegou a 423 MPa no teste). Com a pressão imposta no poço, a pressão na fratura é conhecida, e a abertura pode ser comparada diretamente com Sneddon.

---

## 2. O problema: a abertura derivava com a pressão constante

O caso foi rodado com a pressão no poço mantida em 11 MPa por 6 passos. A solução física não muda nesse intervalo.

![Perfis da abertura](figuras/fig3_perfis.svg)

*Figura 3 — Meia abertura ao longo da fratura. (a) No fim do 1º passo, as duas versões estão próximas de Sneddon. (b) No fim do 6º passo, a versão anterior abriu bem acima da solução, e o excesso cresce em direção à ponta. A versão corrigida não muda.*

![Deriva da abertura no centro](figuras/fig4_deriva.svg)

*Figura 4 — Abertura no centro dividida pela de Sneddon, passo a passo, com a pressão constante.*

| Passo | 1 | 2 | 3 | 4 | 5 | 6 |
|---|---|---|---|---|---|---|
| Antes da correção: erro no centro | +0,35 % | +3,8 % | +6,7 % | +8,8 % | +10,5 % | +11,7 % |
| Corrigido: erro no centro | −1,39 % | −1,39 % | −1,39 % | −1,39 % | −1,39 % | −1,39 % |

No 6º passo, o volume da fratura estava 22 % acima de Sneddon na versão anterior e 2,4 % abaixo na corrigida.

A versão anterior parece melhor no 1º passo (+0,35 %) por acaso. Cada passo tem dois subpassos de acoplamento, e o segundo já soma um pouco da deriva, que compensa a rigidez do domínio finito. Essa compensação deixa de valer a partir do 2º passo.

### 2.1 Causa

Trecho original de `UpdateHydraulicFracture` (`HydraulicFractureModel.f90`), executado a cada subpasso de acoplamento para cada elemento aberto da fratura:

```fortran
sn = <σyy médio nos pontos de integração do elemento> - Biot * p_poros   ! tensão normal total
...
fy = (p + sn) * area * HFmodel%Damage(i)        ! resíduo aplicado como carga na face
IF ( .NOT.HFmodel%FracClosure .AND. fy<0.0 ) fy = 0.0
```

A carga aplicada na face é o **resíduo** entre a pressão do fluido `p` e a tensão normal `sn` do elemento vizinho à face (negativa em compressão). Esse esquema só converge se `sn` tender a −p em todos os elementos da fratura:

- **No primeiro carregamento**, a tensão é a inicial e uniforme (−Shmin), e o resíduo `(p − Shmin)·área` é exatamente a carga de Sneddon.
- **No interior da fratura**, depois de carregado, o elemento fica com σyy ≈ −p, e o resíduo some.
- **Junto à ponta**, a tensão real é singular (cresce como 1/√r). A tensão média de um elemento finito ali não chega a −p, e o resíduo nunca zera. A cada subpasso ele é reaplicado e soma carga nova, então a fratura abre um pouco mais a cada passo.

![Resíduo por elemento](figuras/fig5_residuo.svg)

*Figura 5 — Resíduo `p + σyy médio` em cada elemento da face, no fim do 1º passo (versão anterior). Longe da ponta ele fica abaixo de 33 kPa. Nos três últimos elementos oscila entre −1,9 e +2,4 MPa, e é essa carga espúria que alimenta a deriva.*

---

## 3. A correção

Em `HydraulicFractureModel.f90`, cada elemento da fratura passa a guardar a pressão já aplicada como carga na sua face (`AppliedPressure`) e um indicador de que já foi carregado (`FaceLoaded`):

```fortran
IF( HFmodel%FaceLoaded(i) )THEN
  fy = (p - HFmodel%AppliedPressure(i)) * area * HFmodel%Damage(i)   ! só a variação da pressão
ELSE
  fy = (p + sn) * area * HFmodel%Damage(i)                           ! primeira carga: libera a tensão da face
END IF
IF ( .NOT.HFmodel%FracClosure .AND. fy<0.0 )THEN
  fy = 0.0                                   ! sem fechamento: a queda de pressão não é aplicada
ELSE
  HFmodel%AppliedPressure(i) = p
  HFmodel%FaceLoaded(i) = .TRUE.
END IF
```

Por que funciona:

- **Primeira carga de cada elemento** (quando ele abre, seja pré-fratura ou propagação): continua `p + σyy`. Ela substitui a tração que a rocha transmitia pela face pela pressão do fluido. Para a pré-fratura sob tensão uniforme, é exatamente `(p − Shmin)·área`.
- **Cargas seguintes:** só `Δp·área`. Em elasticidade linear, a superposição das cargas é exata, e o resultado deixa de depender da tensão média do elemento. A singularidade da ponta não realimenta mais a carga.
- **Sem fechamento (`FracClosure` falso):** uma queda de pressão continua não sendo aplicada, como antes, e `AppliedPressure` não é atualizada. Assim uma nova subida de pressão é medida a partir da última carga efetivamente aplicada.
- **Nenhuma carga é contada em dobro:** `UpdateHydraulicFracture` é chamada uma vez por subpasso de acoplamento, antes do laço de iterações da fratura. As re-soluções desse laço reaplicam a mesma carga, e o estado final do subpasso é aceito.

---

## 4. Resultados com a correção

### 4.1 Convergência

O erro em relação a Sneddon tem duas origens: a discretização (elementos lineares, que não representam a raiz quadrada da ponta) e o contorno fixo a distância finita, que deixa o modelo um pouco mais rígido que o meio infinito. As duas fazem o modelo abrir menos que a solução, e as duas diminuem com o refino da malha e com o aumento do domínio:

![Convergência](figuras/fig6_convergencia.svg)

*Figura 6 — Erro relativo (em valor absoluto) em função de dx, para domínios de 10 L e 20 L (versão corrigida, 1º passo). Todos os erros têm sinal negativo.*

| Domínio | dx (m) | Abertura no centro | Perfil (L2) | Volume |
|---|---|---|---|---|
| 10 L (50 m) | 0,5 | −3,29 % | 4,00 % | −5,22 % |
| 10 L (50 m) | 0,25 | −2,45 % | 2,90 % | −3,40 % |
| 10 L (50 m) | 0,125 | −2,02 % | 2,28 % | −2,49 % |
| 20 L (100 m) | 0,5 | −2,26 % | 3,01 % | −4,21 % |
| **20 L (100 m)** | **0,25** | **−1,39 %** | **1,90 %** | **−2,35 %** |
| 20 L (100 m) | 0,125 | −0,94 % | 1,25 % | −1,42 % |
| 40 L (200 m) | 0,25 | −1,12 % | 1,65 % | −2,08 % |

A linha em negrito é a configuração do teste automatizado.

### 4.2 Teste automatizado

[`tests/test_sneddon_injection.py`](../../tests/test_sneddon_injection.py) roda o iGeo e verifica:

| Verificação | Critério | Valor obtido |
|---|---|---|
| Pressão na fratura | uniforme e igual à do poço (±0,5 %) | sem variação nos dígitos gravados; Δp = 1000,02 kPa |
| Abertura no centro | ±2 % de Sneddon | −1,39 % |
| Perfil da abertura | erro L2 < 3 % | 1,90 % |
| Volume da fratura | ±3 % de Sneddon | −2,35 % |
| Face fora da fratura | `uy = 0` | 0 |
| Convergência | o erro diminui de dx = 0,5 m para 0,25 m | sim, nas três métricas |
| Pressão mantida por 3 passos | variação de `uy(0)` < 0,1 % | idêntico em todos os dígitos gravados |

```bash
make -C iGeo_vrs2.0 test               # 7 testes, cerca de 17 s
make -C iGeo_vrs2.0 test MODE=debug    # com verificação de limites, cerca de 70 s
```

A versão anterior falha no teste de pressão mantida. Com o módulo de Young 10 % errado no `.geo`, falham os testes de abertura, perfil e volume, o que confirma que os critérios detectam erros reais.

---

## 5. Efeito em fraturas que se propagam

A correção também muda simulações em que a fratura cresce. O caso [`Simulation_cases/M3_injection_test.geo`](../../Simulation_cases/M3_injection_test.geo) foi rodado com as duas versões. É uma pré-fratura com critério de propagação por dano (`#CDM_FRACTURE`), geomecânica e injeção a vazão crescente, numa malha de 60 mil elementos. O caso simula 360 min e levaria horas, então as simulações foram interrompidas depois de cerca de 70 minutos de processamento.

![Propagação no M3_injection_test](figuras/fig7_propagacao.svg)

*Figura 7 — Pressão máxima no modelo durante o 2º passo de injeção, até 22 min de tempo simulado. Os traços na base marcam os subpassos em que a fratura avançou.*

- No fim do 1º passo (0,1 min), a fratura é idêntica nas duas versões: 1 m de comprimento, topo em −2493,53 m e base em −2495,59 m.
- Até cerca de 2 min, a pressão é a mesma.
- Depois, na **versão anterior**, a pressão sobe de forma suave até 30,5 MPa, e a fratura não volta a se propagar (1 evento, só no início).
- Na **versão corrigida**, a fratura avança 43 vezes em 22 min. A cada avanço a pressão cai e volta a subir (dente de serra), e em 20 min ela está em 27,0 MPa.

A explicação mais provável: na versão anterior, a carga espúria da deriva mantinha a fratura aberta além do equilíbrio elástico, o que aumentava seu volume e absorvia o fluido injetado sem propagar. Como não há solução de referência para este caso, não dá para afirmar qual comportamento está correto só com esta comparação. O resultado mostra que **a correção muda os resultados de propagação** e que eles precisam ser validados.

---

## 6. Outra correção feita no caminho

`ViscoplasticExplicitySolver.f90`: o vetor do pré-condicionador (`diag_precon`) recebia o valor 0 antes de ser alocado, o que é comportamento indefinido no Fortran. Com o compilador Intel isso não fazia nada; com o gfortran causava falha de segmentação no segundo passo estrutural de qualquer caso com fratura. A linha foi removida, porque a rotina seguinte (`CGpreConditioner`) aloca e preenche o vetor.

---

## 7. Limitações e próximos passos

1. **Validar a propagação.** Reproduzir uma fratura KGD (deformação plana, vazão constante, regime de tenacidade ou de viscosidade) e comparar com as soluções de referência (Detournay, 2004; benchmark de Lecampion et al., 2013). O critério de propagação também usa grandezas médias do elemento da ponta (integral J calculada com σyy e εyy do elemento), e pode ter a mesma sensibilidade à malha que foi corrigida aqui.
2. **Vazão controlada em rocha impermeável.** O acoplamento sequencial fluxo → mecânica diverge nesse regime. Uma solução seria incluir a complacência da fratura na equação de fluxo (acoplamento implícito) ou iterar fluxo e mecânica dentro do passo até convergir o volume.
3. **Comparar com os resultados do autor.** Rodar os casos de `Simulation_cases/` com e sem a correção no build Intel original e comparar com os resultados publicados.

---

## 8. Como gerar o PDF e reproduzir as figuras

### 8.1 PDF do relatório

O [`Makefile`](Makefile) desta pasta converte este README em PDF com pandoc e xelatex (MacTeX ou TeX Live):

```bash
# gera docs/relatorio_sneddon/build/relatorio_sneddon.pdf
make -C docs/relatorio_sneddon
```

No PDF, as figuras entram pelas versões `.pdf` de `figuras/`, versionadas junto com as `.svg`, e o filtro [`pdf/figuras_pdf.lua`](pdf/figuras_pdf.lua) junta cada imagem com a sua legenda. As fontes padrão (Palatino e Menlo) são do macOS. Em outro sistema, indique fontes instaladas que tenham os símbolos gregos, por exemplo, no Linux, `make FONTE="DejaVu Serif" FONTE_MONO="DejaVu Sans Mono"`.

### 8.2 Figuras

As figuras são geradas por [`gerar_figuras.py`](gerar_figuras.py) a partir de simulações feitas com o executável atual (`iGeo_vrs2.0/build/iGeo`) e com um executável compilado a partir da versão anterior de `HydraulicFractureModel.f90` (commit `b8828cb`). A docstring do script descreve a estrutura esperada da pasta com as simulações (`<pasta>`):

- o caso do teste com 6 passos de pressão mantida;
- as variantes de malha e de domínio do estudo de convergência;
- o caso `M3_injection_test`.

```bash
make -C docs/relatorio_sneddon figuras DADOS=<pasta>
# equivale a:
.venv/bin/python3 docs/relatorio_sneddon/gerar_figuras.py <pasta>
```

Cada figura é gravada em SVG (usada neste README) e em PDF (usada no PDF do relatório).

---

## Referências

- I. N. Sneddon (1946). The distribution of stress in the neighbourhood of a crack in an elastic solid. *Proc. R. Soc. Lond. A* 187, 229–260.
- I. N. Sneddon, M. Lowengrub (1969). *Crack Problems in the Classical Theory of Elasticity*. Wiley.
- E. Detournay (2004). Propagation regimes of fluid-driven fractures in impermeable rocks. *Int. J. Geomech.* 4(1), 35–45.
- B. Lecampion et al. (2013). The impact of the near-tip logic on the accuracy and convergence rate of hydraulic fracture simulators compared to reference solutions. Em *Effective and Sustainable Hydraulic Fracturing*, InTech. [PDF](https://personal.math.ubc.ca/~peirce/Lecampion_et_al_BenchmarkPaper_2013.pdf)
- E. R. Simonson, A. S. Abou-Sayed, R. J. Clifton (1978). Containment of massive hydraulic fractures. *SPE Journal* 18(1), 27–32. [OSTI](https://www.osti.gov/biblio/5769471)
