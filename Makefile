PYTHON ?= python3

.PHONY: install install-dev metrics cluster test reproduce

install:
	$(PYTHON) -m pip install -r requirements.txt

install-dev:
	$(PYTHON) -m pip install -r requirements-dev.txt

metrics:
	$(PYTHON) -m grav.metrics --out-dir results

cluster:
	$(PYTHON) -m grav.cluster --csv results/eca_metrics.csv --out-dir results

test:
	$(PYTHON) -m unittest discover -s tests -p 'test_*.py'

reproduce: metrics cluster
