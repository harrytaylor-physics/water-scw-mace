"""Root pytest configuration: keep `pytest` from the repository root collecting
this repository's suite (tests/) only, as before the subneptune migration.

- production_campaign/ holds the ARC scripts verbatim, as provenance; one of
  them (scw/residual_lag_test.py) matches pytest's *_test.py pattern and parses
  sys.argv at import, so collecting it aborts the run.
- subneptune/ has its own suite and configuration; run it from there
  (`cd subneptune && python -m pytest`).
"""

collect_ignore = ["production_campaign", "subneptune"]
