"""Count saved repetition results without loading models or requiring R packages."""
import argparse
from datetime import datetime
from pathlib import Path
import time

METHODS = ("SurvDeepFIT", "Cox", "AFT", "RSF", "XGBoost")


def report(runs, repetitions):
    total_finished = total_started = total_records = 0
    rows = []
    for scenario in (1, 2):
        for n in (1000, 3000):
            for rho in ("0", "0.3", "0.5"):
                setting = runs / f"scenario{scenario}" / f"n{n}" / f"rho{rho}"
                finished = started = records = 0
                for directory in setting.glob("rep*"):
                    if not directory.name[3:].isdigit():
                        continue
                    if not 1 <= int(directory.name[3:]) <= repetitions:
                        continue
                    started += (directory / "run.rds").is_file()
                    saved = sum((directory / method / "result.rds").is_file()
                                for method in METHODS)
                    records += saved
                    finished += saved == len(METHODS)
                rows.append(f"{scenario:>8} {n:>6} {rho:>5} {finished:>9}/{repetitions:<5}")
                total_finished += finished
                total_started += started
                total_records += records
    planned = 12 * repetitions
    print(f"\n{datetime.now().astimezone().isoformat(timespec='seconds')}")
    print(f"Results: {runs}")
    print(f"Finished repetitions: {total_finished:,}/{planned:,} "
          f"({100 * total_finished / planned:.1f}%)")
    print(f"Started (run metadata saved): {total_started:,}/{planned:,}")
    print(f"Saved method results: {total_records:,}/{planned * len(METHODS):,}")
    print("Scenario      n   rho  Finished / planned")
    print("\n".join(rows))
    print("Finished = all five result.rds files saved.")
    print("Saved results can include failed fits; this counts completion, not fit success.",
          flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path, help="Study directory or its runs directory")
    parser.add_argument("--repetitions", type=int, default=1000,
                        help="Planned repetitions per setting (default: 1000)")
    parser.add_argument("--watch", type=float, metavar="SECONDS",
                        help="Refresh interval; Ctrl+C stops monitoring, leaving jobs running")
    args = parser.parse_args()
    if args.repetitions < 1 or (args.watch is not None and args.watch < 10):
        parser.error("Use positive repetitions and a watch interval of at least 10 seconds.")
    directory = args.directory.resolve()
    if not directory.is_dir():
        parser.error(f"Directory does not exist: {directory}")
    runs = directory / "runs" if (directory / "runs").is_dir() else directory
    try:
        while True:
            report(runs, args.repetitions)
            if args.watch is None:
                break
            time.sleep(args.watch)
    except KeyboardInterrupt:
        pass


if __name__ == "__main__":
    main()
