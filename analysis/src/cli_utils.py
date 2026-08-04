"""
cli_utils.py
============
Shared argparse/logging/timing boilerplate for analysis pipeline scripts.

Each script builds its own ArgumentParser and adds its own flags; this only
consolidates what is identical across all of them.
"""

import argparse
import os
import time

from utils import setup_logging


def add_common_args(parser: argparse.ArgumentParser, *, force: bool = False) -> None:
    """Register --verbose, and --force where the script supports it."""
    parser.add_argument('--verbose', action='store_true',
                        help='Show INFO-level messages on console')
    if force:
        parser.add_argument('--force', action='store_true',
                            help='Re-run even if output already exists')


def init_run(module_name: str, args: argparse.Namespace, log_dir: str) -> float:
    """log setup + start banner; returns t0 for finish_run()."""
    log_path = os.path.join(log_dir, f'{module_name}.log')
    setup_logging(args.verbose, log_path)
    print(f'=== {module_name} ===')
    return time.time()


def finish_run(t0: float) -> None:
    print(f'\nDone in {time.time()-t0:.1f}s.')
