# Reproduce the manuscript from the committed data.
#
# `make all` runs the chain in dependency order: derived tables, then the
# register of every number the paper quotes, then the figures, then the two
# audits that check the register against the manuscript in both directions.

PY := python3
PAPER := paper

.PHONY: all numbers stage1 stage2 figures audit test paper clean help

help:
	@echo "make numbers   recompute the derived tables and the number register"
	@echo "make figures   redraw Figures 1-5 from the derived tables"
	@echo "make audit     check every number in the register against paper.tex"
	@echo "make test      run the test suite"
	@echo "make paper     build paper.pdf (needs pdflatex and bibtex)"
	@echo "make all       numbers, figures, audit, test"
	@echo "make clean     remove build artefacts, keep committed outputs"

all: numbers figures audit test

# ---------------------------------------------------------------- numbers --
# Each script writes one or more tables under data/derived/. The register is
# built last because it reads all of them.
numbers: stage1 stage2

# Stage 1 builds the per-gene decomposition tables the rest of the pipeline
# reads. These live under notebook/ for historical reasons: they are the
# original analysis, kept as the record of how the decompositions were first
# derived. They are plain scripts, not notebooks, and run in order.
stage1:
	$(PY) notebook/01_t2_supply_table.py
	$(PY) notebook/03_clinical_cohorts.py
	$(PY) notebook/04_oh_decomposition.py
	$(PY) notebook/05_prospective_predictions.py
	$(PY) notebook/06_compound_predictions.py
	$(PY) notebook/07_kras_g12c_prospective.py

# Stage 2 is everything downstream of the decompositions. The register is
# built last because it reads all of them.
stage2:
	$(PY) scripts/recompute_topk_with_ci.py
	$(PY) scripts/recompute_benchmark_r3.py
	$(PY) scripts/recompute_enu_class_test.py
	$(PY) scripts/recompute_retrospective_validation.py
	$(PY) scripts/two_hit_reachability.py
	$(PY) scripts/sigma_transport_validation.py
	$(PY) scripts/oh_dose_transport.py
	$(PY) scripts/met_audit.py
	$(PY) scripts/build_reachability_track.py
	$(PY) scripts/env_versions.py
	$(PY) scripts/compute_paper_numbers.py
	$(PY) scripts/build_pipelines_md.py
# Last, so the checksums fix the tables as this run left them. Run by hand it
# goes stale the moment anything upstream is recomputed.
	$(PY) scripts/sha256sum.py

# ---------------------------------------------------------------- figures --
figures:
	$(PY) $(PAPER)/figure_1.py
	$(PY) $(PAPER)/figure_2.py
	$(PY) $(PAPER)/figure_3.py
	$(PY) $(PAPER)/figure_4.py
	$(PY) $(PAPER)/figure_5.py

# ------------------------------------------------------------------ audit --
# Two directions. The first checks that every registered value appears in the
# manuscript where the register says it does; the second checks that every
# number written in the prose is accounted for by a registered value. The
# second is the one that catches a hand-typed figure going stale.
audit:
	$(PY) scripts/audit_paper_numbers.py
	$(PY) scripts/audit_unregistered_numbers.py

test:
	$(PY) -m pytest tests/ -q

# ------------------------------------------------------------------ paper --
paper:
	cd $(PAPER) && pdflatex -interaction=nonstopmode paper.tex >/dev/null
	cd $(PAPER) && bibtex paper >/dev/null
	cd $(PAPER) && pdflatex -interaction=nonstopmode paper.tex >/dev/null
	cd $(PAPER) && pdflatex -interaction=nonstopmode paper.tex >/dev/null
	@echo "wrote $(PAPER)/paper.pdf"

clean:
	find . -name '__pycache__' -type d -exec rm -rf {} + 2>/dev/null || true
	rm -rf .pytest_cache
	rm -f $(PAPER)/*.aux $(PAPER)/*.log $(PAPER)/*.out \
	      $(PAPER)/*.bbl $(PAPER)/*.blg
