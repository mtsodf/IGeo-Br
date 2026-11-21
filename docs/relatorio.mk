# Regras comuns para gerar os relatórios de docs/ em PDF.
#
# Cada relatório tem um Makefile que define NOME (nome do PDF) e inclui este
# arquivo. Requer pandoc e um LaTeX com xelatex (MacTeX/TeX Live). As figuras
# entram no PDF pelas versões .pdf de figuras/, versionadas junto com as .svg;
# o filtro e os metadados do pandoc ficam em docs/pdf/.

DOCS       := $(patsubst %/,%,$(dir $(abspath $(lastword $(MAKEFILE_LIST)))))
PANDOC     ?= pandoc
PDF_ENGINE ?= xelatex
# Ambiente da raiz (requirements.txt) se existir; senão, python3 do PATH.
VENV_PYTHON := $(DOCS)/../.venv/bin/python3
PYTHON     ?= $(if $(wildcard $(VENV_PYTHON)),$(abspath $(VENV_PYTHON)),python3)
FONTE      ?= Palatino
FONTE_MONO ?= Menlo
DATA       ?= $(shell date +%d/%m/%Y)

BUILD    := build
SAIDA    := $(BUILD)/$(NOME).pdf
FONTE_MD := README.md
FIGURAS  := $(patsubst %.svg,%.pdf,$(wildcard figuras/*.svg))
APOIO    := $(DOCS)/pdf/metadados.yaml $(DOCS)/pdf/figuras_pdf.lua $(DOCS)/relatorio.mk Makefile

.PHONY: all pdf figuras clean help

all: pdf

pdf: $(SAIDA)

# --shift-heading-level-by=-1: o título "# ..." do README vira o título do PDF.
# -implicit_figures: a legenda é o parágrafo em itálico que segue cada imagem.
$(SAIDA): $(FONTE_MD) $(FIGURAS) $(APOIO) | $(BUILD)
	$(PANDOC) $(FONTE_MD) \
	  --from=markdown-implicit_figures \
	  --shift-heading-level-by=-1 \
	  --metadata-file=$(DOCS)/pdf/metadados.yaml \
	  --metadata=date:"$(DATA)" \
	  --variable=mainfont:"$(FONTE)" \
	  --variable=monofont:"$(FONTE_MONO)" \
	  --variable=monofontoptions:Scale=0.85 \
	  --lua-filter=$(DOCS)/pdf/figuras_pdf.lua \
	  --resource-path=. \
	  --pdf-engine=$(PDF_ENGINE) \
	  -o $@
	@echo "PDF gerado: $(CURDIR)/$@"

$(BUILD):
	mkdir -p $@

# Sem pré-requisitos: só roda se a figura PDF não existir.
figuras/%.pdf:
	@echo "Falta $@. Regenere as figuras com: make figuras DADOS=<pasta_com_as_simulações>" >&2
	@exit 1

figuras:
	@test -n "$(DADOS)" || { echo "Uso: make figuras DADOS=<pasta_com_as_simulações>" >&2; exit 1; }
	$(PYTHON) gerar_figuras.py "$(DADOS)"

clean:
	rm -f $(SAIDA)

help:
	@sed -n '2,8p' Makefile
