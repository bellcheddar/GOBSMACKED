# Small, deliberately: everything else is a pixi or uv task.
PY := .venv/bin/python

.PHONY: test check-refs third-party serve

test:
	$(PY) -m pytest tests -q

# Every DOI in software.yaml must resolve. A reference that has not been checked
# does not ship.
check-refs:
	$(PY) scripts/check_refs.py

# THIRD_PARTY.md is generated, never hand-edited: software.yaml is the source.
third-party:
	$(PY) scripts/gen_third_party.py

serve:
	$(PY) -m flask --app app run --debug --port 8009
