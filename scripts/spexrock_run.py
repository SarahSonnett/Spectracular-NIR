"""Source-tree shim for ``spexrock-run`` (no install needed)."""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from spexrock.cli.run import main  # noqa: E402

if __name__ == "__main__":
    sys.exit(main())
