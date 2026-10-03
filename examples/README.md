# Example configs

- `uspex_lxd_sampledata.json` — reduces the pyspextool sample LXD data set
  (a star pair; clone `github.com/pyspextool/test_data` into
  `data/sample/test_data` first). Run from the repo root:

  ```bash
  .venv/bin/python scripts/spexrock_run.py examples/uspex_lxd_sampledata.json
  ```

- `asteroid_template.json` — a starting point for a real asteroid + solar
  analog reduction: fill in paths, frame numbers, and (for LXD thermal
  work) the `neatm_*` circumstances from Horizons.
