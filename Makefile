# LiftOn dev targets. Override the interpreter if your env lives elsewhere:
#   make benchmark-gate LIFTON_PY=/path/to/envs/lifton_devel/bin/python
LIFTON_PY ?= /home/kh.chao/miniconda3/envs/lifton_devel/bin/python

.PHONY: test test-fast benchmark-gate benchmark-gate-update benchmark-suite help

help:
	@echo "LiftOn make targets:"
	@echo "  test                  run the full pytest suite"
	@echo "  test-fast             run the 24-cell byte-identity gate + integration tests"
	@echo "  benchmark-gate        pytest gate + fast single-chrom benchmark vs baseline"
	@echo "  benchmark-gate-update reseed the fast-benchmark baseline (after a reviewed change)"
	@echo "  benchmark-suite       run the full 8-benchmark comparison + report (on demand)"

# Full suite. The 3 hypothesis-only files error on collection without hypothesis
# installed; deselect them for a clean run (see CLAUDE.md).
test:
	PYTHONNOUSERSITE=1 $(LIFTON_PY) -m pytest tests/ -q \
		--ignore=tests/test_property_based.py \
		--ignore=tests/test_streaming_property.py \
		--ignore=tests/test_vulnerabilities.py

test-fast:
	PYTHONNOUSERSITE=1 $(LIFTON_PY) -m pytest \
		tests/test_native_matrix.py tests/test_integration_pipeline.py -q

benchmark-gate:
	PYTHONNOUSERSITE=1 $(LIFTON_PY) scripts/benchmark_gate.py

benchmark-gate-update:
	PYTHONNOUSERSITE=1 $(LIFTON_PY) scripts/benchmark_gate.py --update-baseline

benchmark-suite:
	PYTHONNOUSERSITE=1 $(LIFTON_PY) -m benchmarks.compare.run_compare --all -t 8 -j 2
