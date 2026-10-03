"""``spexrock-run``: reduce one SpeX data set from a JSON config.

Typical use::

    spexrock-run --init my_asteroid.json   # write a template config
    <edit my_asteroid.json>
    spexrock-run my_asteroid.json          # run the pipeline
    spexrock-run my_asteroid.json --fresh  # ignore existing intermediates

The template config carries every field of
:class:`~spexrock.config.ReduceConfig` with its default value, so editing it
is the whole configuration step -- there are no other knobs.
"""

from __future__ import annotations

import argparse
import os
import sys

os.environ.setdefault("MPLBACKEND", "Agg")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="spexrock-run",
        description="Asteroid reflectance reduction for IRTF SpeX "
                    "(pyspextool engine).")
    parser.add_argument("config", help="JSON config file (see --init)")
    parser.add_argument("--init", action="store_true",
                        help="write a template config to CONFIG and exit")
    parser.add_argument("--fresh", action="store_true",
                        help="rerun all stages even if outputs exist")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)

    from spexrock.config import ReduceConfig

    if args.init:
        if os.path.exists(args.config):
            print(f"Refusing to overwrite existing file {args.config}")
            return 1
        ReduceConfig().save(args.config)
        print(f"Wrote template config to {args.config}; edit it, then run\n"
              f"  spexrock-run {args.config}")
        return 0

    config = ReduceConfig.load(args.config)

    from spexrock import pipeline
    pipeline.run(config, resume=not args.fresh)
    return 0


if __name__ == "__main__":
    sys.exit(main())
