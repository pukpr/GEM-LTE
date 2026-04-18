#!/usr/bin/env python3

from __future__ import annotations

import argparse
import errno
import json
import os
import pty
import re
import select
import shlex
import shutil
import subprocess
import termios
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any


DEFAULT_RESULTS_HTML = Path(__file__).resolve().parents[1] / "docs" / "gem-lte-results.html"
STATUS_RE = re.compile(
    r"Status:\s+"
    r"(?P<thread>\d+)\s+"
    r"(?P<counter>\d+)\s+"
    r"(?P<interval>[-+]?\d+(?:\.\d*)?(?:[Ee][-+]?\d+)?)\s+"
    r"(?P<complement>[-+]?\d+(?:\.\d*)?(?:[Ee][-+]?\d+)?)"
    r"(?:\s+#\s*(?P<loop>\d+))?\s*$"
)
FINAL_RE = re.compile(
    r"(?P<label>[A-Za-z][A-Za-z0-9_]*)\s+"
    r"(?P<interval>[-+]?\d+(?:\.\d*)?(?:[Ee][-+]?\d+)?)\s+"
    r"(?P<complement>[-+]?\d+(?:\.\d*)?(?:[Ee][-+]?\d+)?)\s+"
    r"(?P<thread>\d+)\s+"
    r"(?P<counter>\d+)\s*$"
)


@dataclass
class MetricEvent:
    kind: str
    interval_cc: float
    complement_cc: float
    thread: int
    counter: int
    loop_index: int | None = None
    label: str | None = None
    line: str | None = None


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Run lt.exe/lt over Feb2026 experiment directories, watch the live "
            "Status stream, and inject keys when CC thresholds are reached."
        )
    )
    parser.add_argument("experiments_root", type=Path, help="Path such as experiments/Feb2026")
    parser.add_argument(
        "evaluation_root",
        nargs="?",
        type=Path,
        help=(
            "Mirror workspace where copied run directories, logs, and JSON summaries are written "
            "(default: sibling of experiments_root named April2026-expect)"
        ),
    )
    parser.add_argument(
        "max_runs",
        nargs="?",
        type=int,
        default=None,
        help="Maximum number of ranked runs to process",
    )
    parser.add_argument("--results-html", type=Path, default=DEFAULT_RESULTS_HTML)
    parser.add_argument("--site", action="append", dest="sites", default=[], help="Run only this site ID")
    parser.add_argument("--threshold", type=float, default=0.6, help="Threshold to watch for (default: 0.6)")
    parser.add_argument(
        "--watch",
        choices=["test", "train", "interval", "complement"],
        default="test",
        help="Which live metric to watch (default: test)",
    )
    parser.add_argument(
        "--action",
        choices=["save", "quit", "digit"],
        default="save",
        help=(
            "Key to send once the watched metric crosses the threshold "
            "(default: save). With digit, lt.exe uses the enso_opt mapping "
            "1->0.1, 2->0.2, ..., 9->0.9."
        ),
    )
    parser.add_argument(
        "--trigger-key",
        help="Initial digit to send when --action=digit. Defaults to the threshold tenths digit (for example 0.6 -> 6).",
    )
    parser.add_argument(
        "--trigger-floor",
        type=float,
        help=(
            "Lowest digit-trigger threshold allowed for --action=digit. "
            "Defaults to two tenths below --threshold when possible (for example 0.6 -> 0.4)."
        ),
    )
    parser.add_argument(
        "--dynamic-trigger",
        action=argparse.BooleanOptionalAction,
        default=True,
        help=(
            "With --action=digit, set the initial trigger immediately and lower it "
            "toward --trigger-floor as TIMEOUT elapses if the watched CC still "
            "has not reached --threshold."
        ),
    )
    parser.add_argument("--metric", default="CC")
    parser.add_argument("--timeout", type=float, default=1000.0)
    parser.add_argument("--train-start", default="1940")
    parser.add_argument("--train-end", default="1970")
    parser.add_argument("--idate", default="1920.9")
    parser.add_argument("--json", dest="use_json", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--exclude", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--filter", dest="filter9", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--trend", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--test-only", action=argparse.BooleanOptionalAction, default=False)
    parser.add_argument("--sim", action=argparse.BooleanOptionalAction, default=False)
    parser.add_argument("--reset", action="store_true", help="Delete evaluation_root before running")
    return parser.parse_args()


def parse_ranked_sites(results_html: Path) -> list[dict[str, Any]]:
    html = results_html.read_text(encoding="utf-8")
    match = re.search(r"const ALL_SITES = (\[[\s\S]*?\]);", html)
    if not match:
        raise RuntimeError(f"ALL_SITES not found in {results_html}")
    sites = json.loads(match.group(1))
    sites.sort(key=lambda site: -float(site["r_val"]))
    return sites


def default_evaluation_root(experiments_root: Path) -> Path:
    return (experiments_root.parent / "April2026-expect").resolve()


def resolve_lt_cmd(root: Path, use_json: bool) -> list[str]:
    for name in ("lt.exe", "lt"):
        candidate = root / name
        if candidate.exists():
            cmd = [str(candidate)]
            if use_json:
                cmd.append("-j")
            return cmd
    raise FileNotFoundError(f"No lt.exe or lt found in {root}")


def available_sites(args: argparse.Namespace) -> list[dict[str, Any]]:
    if args.sites:
        selected = []
        ranked = {str(site["id"]): site for site in parse_ranked_sites(args.results_html)} if args.results_html.exists() else {}
        for site_id in args.sites:
            site_dir = args.experiments_root / site_id
            if not site_dir.is_dir():
                raise FileNotFoundError(f"Site directory not found: {site_dir}")
            selected.append(ranked.get(site_id, {"id": site_id, "name": site_id, "r_val": 0.0}))
        return selected

    if args.results_html.exists():
        ranked = parse_ranked_sites(args.results_html)
        usable = []
        for site in ranked:
            site_id = str(site["id"])
            site_dir = args.experiments_root / site_id
            if site_dir.is_dir() and (site_dir / "lt.exe.p").is_file() and (site_dir / f"{site_id}.dat").is_file():
                usable.append(site)
        return usable

    usable = []
    for site_dir in sorted(args.experiments_root.iterdir()):
        if not site_dir.is_dir():
            continue
        site_id = site_dir.name
        if (site_dir / "lt.exe.p").is_file() and (site_dir / f"{site_id}.dat").is_file():
            usable.append({"id": site_id, "name": site_id, "r_val": 0.0})
    return usable


def metric_names(exclude: bool) -> tuple[str, str]:
    if exclude:
        return "train", "test"
    return "test", "train"


def event_semantics(event: MetricEvent, exclude: bool) -> dict[str, float]:
    primary_name, secondary_name = metric_names(exclude)
    return {
        "interval_cc": event.interval_cc,
        "complement_cc": event.complement_cc,
        f"{primary_name}_cc": event.interval_cc,
        f"{secondary_name}_cc": event.complement_cc,
    }


def watched_value(event: MetricEvent, watch: str, exclude: bool) -> float:
    if watch == "interval":
        return event.interval_cc
    if watch == "complement":
        return event.complement_cc
    primary_name, secondary_name = metric_names(exclude)
    if watch == primary_name:
        return event.interval_cc
    if watch == secondary_name:
        return event.complement_cc
    raise ValueError(f"Unsupported watch target: {watch}")


def infer_trigger_key(threshold: float) -> str:
    scaled = round(threshold * 10.0)
    if abs(scaled / 10.0 - threshold) > 1.0e-9 or not (1 <= scaled <= 9):
        raise ValueError(
            "--action=digit requires --trigger-key or a threshold aligned with digits 1..9 "
            "(exactly as enso_opt maps 6 -> 0.6, 7 -> 0.7, and so on)"
        )
    return str(scaled)


def infer_trigger_digit(value: float, option_name: str) -> int:
    scaled = round(value * 10.0)
    if abs(scaled / 10.0 - value) > 1.0e-9 or not (1 <= scaled <= 9):
        raise ValueError(f"{option_name} must align to one of the lt.exe trigger digits 0.1..0.9")
    return scaled


def resolve_digit_trigger_policy(args: argparse.Namespace) -> tuple[int, int]:
    if args.trigger_key is not None:
        if len(args.trigger_key) != 1 or args.trigger_key < "1" or args.trigger_key > "9":
            raise ValueError("--trigger-key must be a single digit 1..9")
        initial_digit = int(args.trigger_key)
    else:
        initial_digit = infer_trigger_digit(args.threshold, "--threshold")

    floor_digit = (
        infer_trigger_digit(args.trigger_floor, "--trigger-floor")
        if args.trigger_floor is not None
        else max(1, initial_digit - 2)
    )
    if floor_digit > initial_digit:
        raise ValueError("--trigger-floor must be less than or equal to the initial digit trigger")
    return initial_digit, floor_digit


def parse_event(line: str) -> MetricEvent | None:
    stripped = line.strip()
    if not stripped:
        return None
    status_match = STATUS_RE.match(stripped)
    if status_match:
        return MetricEvent(
            kind="status",
            interval_cc=float(status_match.group("interval")),
            complement_cc=float(status_match.group("complement")),
            thread=int(status_match.group("thread")),
            counter=int(status_match.group("counter")),
            loop_index=int(status_match.group("loop")) if status_match.group("loop") else None,
            line=stripped,
        )
    final_match = FINAL_RE.match(stripped)
    if final_match:
        return MetricEvent(
            kind="final",
            interval_cc=float(final_match.group("interval")),
            complement_cc=float(final_match.group("complement")),
            thread=int(final_match.group("thread")),
            counter=int(final_match.group("counter")),
            label=final_match.group("label"),
            line=stripped,
        )
    return None


def configure_slave_tty(fd: int) -> None:
    attrs = termios.tcgetattr(fd)
    attrs[3] &= ~(termios.ECHO | termios.ICANON)
    attrs[6][termios.VMIN] = 1
    attrs[6][termios.VTIME] = 0
    termios.tcsetattr(fd, termios.TCSANOW, attrs)


def send_key(master_fd: int, key: str) -> bool:
    try:
        os.write(master_fd, key.encode("ascii"))
        return True
    except OSError as exc:
        if exc.errno == errno.EIO:
            return False
        raise


def create_sim_data(run_dir: Path) -> str:
    csv_path = run_dir / "lte_results.csv"
    if not csv_path.exists():
        raise FileNotFoundError(f"lte_results.csv not found in {run_dir}")
    sim_dat = run_dir / "lte_results_model.dat"
    with csv_path.open(encoding="utf-8") as src, sim_dat.open("w", encoding="utf-8") as dst:
        for raw_line in src:
            line = raw_line.strip()
            if not line:
                continue
            parts = [part.strip() for part in line.split(",")]
            if len(parts) >= 2:
                dst.write(f"{parts[0]}  {parts[1]}\n")
    return sim_dat.name


def build_env(args: argparse.Namespace, site_id: str, run_dir: Path) -> dict[str, str]:
    env = os.environ.copy()
    env["METRIC"] = args.metric.strip().upper()
    env["TIMEOUT"] = f"{args.timeout:.6f}"
    env["TRAIN_START"] = args.train_start
    env["TRAIN_END"] = args.train_end
    env["CLIMATE_INDEX"] = create_sim_data(run_dir) if args.sim else f"{site_id}.dat"
    env["IDATE"] = args.idate
    env["EXCLUDE"] = str(args.exclude).lower()
    env["TREND"] = str(args.trend).lower()
    env["F9"] = "1" if args.filter9 else "0"
    env["TEST_ONLY"] = str(args.test_only).lower()
    return env


def max_or_none(values: list[float]) -> float | None:
    return max(values) if values else None


def iso_now() -> str:
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())


def mirror_runtime_root(source_root: Path, evaluation_root: Path) -> None:
    evaluation_root.mkdir(parents=True, exist_ok=True)
    for entry in source_root.iterdir():
        if entry.is_file():
            shutil.copy2(entry, evaluation_root / entry.name)
        elif entry.is_dir() and entry.name == "locs":
            target_dir = evaluation_root / entry.name
            if target_dir.exists():
                shutil.rmtree(target_dir)
            shutil.copytree(entry, target_dir)


def mirror_site_dir(source_root: Path, evaluation_root: Path, site_id: str) -> Path:
    source_dir = (source_root / site_id).resolve()
    if not source_dir.is_dir():
        raise FileNotFoundError(f"Site directory not found: {source_dir}")
    target_dir = (evaluation_root / site_id).resolve()
    if target_dir.exists():
        shutil.rmtree(target_dir)
    shutil.copytree(source_dir, target_dir, dirs_exist_ok=True)
    return target_dir


def run_site(
    args: argparse.Namespace,
    lt_cmd: list[str],
    site: dict[str, Any],
    source_root: Path,
    evaluation_root: Path,
) -> dict[str, Any]:
    site_id = str(site["id"])
    site_name = str(site.get("name", site_id))
    source_dir = (source_root / site_id).resolve()
    run_dir = mirror_site_dir(source_root, evaluation_root, site_id)
    log_dir = evaluation_root / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    log_file = log_dir / f"{site_id}.log"
    env = build_env(args, site_id, run_dir)

    par_path = run_dir / f"lt.exe.{site_id}.dat.par"
    if par_path.exists():
        par_path.unlink()

    primary_name, secondary_name = metric_names(args.exclude)
    initial_trigger_digit = None
    trigger_floor_digit = None
    if args.action == "digit":
        initial_trigger_digit, trigger_floor_digit = resolve_digit_trigger_policy(args)
        action_key = str(initial_trigger_digit)
    else:
        action_key = {"save": "s", "quit": "q"}[args.action]
    action_note = None
    if args.action == "digit":
        action_note = (
            f"Digit triggers lt.exe's internal stop rule on the complement metric ({secondary_name}). "
            f"Initial trigger={initial_trigger_digit / 10.0:.1f}, "
            f"floor={trigger_floor_digit / 10.0:.1f}, "
            f"dynamic={'true' if args.dynamic_trigger else 'false'}."
        )

    command = "ulimit -s unlimited && exec " + " ".join(shlex.quote(part) for part in lt_cmd)
    master_fd, slave_fd = pty.openpty()
    configure_slave_tty(slave_fd)
    proc = subprocess.Popen(
        ["bash", "-lc", command],
        cwd=str(run_dir),
        env=env,
        stdin=slave_fd,
        stdout=slave_fd,
        stderr=slave_fd,
        close_fds=True,
    )
    os.close(slave_fd)

    started_at = iso_now()
    started_monotonic = time.monotonic()
    best_events: list[MetricEvent] = []
    final_event: MetricEvent | None = None
    threshold_event: MetricEvent | None = None
    action_sent = False
    best_watched = float("-inf")
    trigger_history: list[dict[str, Any]] = []
    current_trigger_digit = initial_trigger_digit
    buffer = ""
    threshold_notice = f"[threshold] site={site_id} watch={args.watch} value={{value:.6f}} threshold={args.threshold:.6f}"

    with log_file.open("w", encoding="utf-8") as log:
        try:
            if args.action == "digit":
                if send_key(master_fd, action_key):
                    notice = (
                        f"[trigger-set] site={site_id} key={action_key} "
                        f"cc={current_trigger_digit / 10.0:.1f}"
                    )
                    trigger_history.append(
                        {
                            "elapsed_seconds": 0.0,
                            "key": action_key,
                            "trigger_cc": current_trigger_digit / 10.0,
                            "reason": "initial",
                        }
                    )
                    print(notice, flush=True)
                    log.write(notice + "\n")
                    log.flush()
            while True:
                ready, _, _ = select.select([master_fd], [], [], 0.25)
                if (
                    args.action == "digit"
                    and args.dynamic_trigger
                    and current_trigger_digit is not None
                    and trigger_floor_digit is not None
                    and current_trigger_digit > trigger_floor_digit
                    and threshold_event is None
                    and args.timeout > 0.0
                ):
                    elapsed = time.monotonic() - started_monotonic
                    span = current_trigger_digit - trigger_floor_digit
                    total_span = initial_trigger_digit - trigger_floor_digit
                    if total_span > 0:
                        desired_drop = int(min(max(elapsed / args.timeout, 0.0), 1.0) * total_span)
                        desired_digit = max(trigger_floor_digit, initial_trigger_digit - desired_drop)
                        if desired_digit < current_trigger_digit:
                            next_key = str(desired_digit)
                            if send_key(master_fd, next_key):
                                current_trigger_digit = desired_digit
                                notice = (
                                    f"[trigger-lower] site={site_id} key={next_key} "
                                    f"cc={current_trigger_digit / 10.0:.1f} "
                                    f"elapsed={elapsed:.1f}s best_watch="
                                    f"{best_watched if best_watched > float('-inf') else float('nan'):.6f}"
                                )
                                trigger_history.append(
                                    {
                                        "elapsed_seconds": elapsed,
                                        "key": next_key,
                                        "trigger_cc": current_trigger_digit / 10.0,
                                        "reason": "timeout-progress",
                                    }
                                )
                                print(notice, flush=True)
                                log.write(notice + "\n")
                                log.flush()
                if ready:
                    try:
                        chunk = os.read(master_fd, 4096)
                    except OSError as exc:
                        if exc.errno == errno.EIO:
                            if proc.poll() is not None:
                                break
                            continue
                        raise
                    if not chunk:
                        if proc.poll() is not None:
                            break
                        continue
                    text = chunk.decode("utf-8", errors="replace")
                    log.write(text)
                    log.flush()
                    buffer += text.replace("\r\n", "\n").replace("\r", "\n")
                    while "\n" in buffer:
                        line, buffer = buffer.split("\n", 1)
                        event = parse_event(line)
                        if not event:
                            continue
                        best_events.append(event)
                        if event.kind == "final":
                            final_event = event
                        watched = watched_value(event, args.watch, args.exclude)
                        best_watched = max(best_watched, watched)
                        if threshold_event is None and watched >= args.threshold:
                            threshold_event = event
                            notice = threshold_notice.format(value=watched)
                            print(notice, flush=True)
                            log.write(notice + "\n")
                            log.flush()
                            if args.action != "digit" and not action_sent:
                                action_sent = send_key(master_fd, action_key)

                if proc.poll() is not None and not ready:
                    break
        finally:
            try:
                os.close(master_fd)
            except OSError:
                pass

    exit_code = proc.wait()
    interval_values = [event.interval_cc for event in best_events]
    complement_values = [event.complement_cc for event in best_events]
    best_interval = max_or_none(interval_values)
    best_complement = max_or_none(complement_values)

    summary = {
        "site_id": site_id,
        "site_name": site_name,
        "seed_r_val": float(site.get("r_val", 0.0)),
        "source_dir": str(source_dir),
        "run_dir": str(run_dir),
        "log_file": str(log_file.resolve()),
        "started_at": started_at,
        "completed_at": iso_now(),
        "exit_code": exit_code,
        "watch": args.watch,
        "threshold": args.threshold,
        "action": args.action,
        "action_key": action_key,
        "action_note": action_note,
        "dynamic_trigger": args.dynamic_trigger if args.action == "digit" else None,
        "trigger_floor": None if trigger_floor_digit is None else trigger_floor_digit / 10.0,
        "trigger_history": trigger_history,
        "env": {
            "METRIC": env["METRIC"],
            "TIMEOUT": env["TIMEOUT"],
            "TRAIN_START": env["TRAIN_START"],
            "TRAIN_END": env["TRAIN_END"],
            "CLIMATE_INDEX": env["CLIMATE_INDEX"],
            "IDATE": env["IDATE"],
            "EXCLUDE": env["EXCLUDE"],
            "TREND": env["TREND"],
            "F9": env["F9"],
            "TEST_ONLY": env["TEST_ONLY"],
        },
        "metric_roles": {
            "interval": primary_name,
            "complement": secondary_name,
        },
        "trigger_mapping": "digit N maps to cc=N/10 via enso_opt -> Set_Trigger",
        "best_interval_cc": best_interval,
        "best_complement_cc": best_complement,
        f"best_{primary_name}_cc": best_interval,
        f"best_{secondary_name}_cc": best_complement,
        "threshold_crossed": threshold_event is not None,
        "best_watched_cc": None if best_watched == float("-inf") else best_watched,
        "threshold_event": None if threshold_event is None else {
            **asdict(threshold_event),
            **event_semantics(threshold_event, args.exclude),
            "watched_cc": watched_value(threshold_event, args.watch, args.exclude),
        },
        "final_event": None if final_event is None else {
            **asdict(final_event),
            **event_semantics(final_event, args.exclude),
        },
        "event_count": len(best_events),
    }
    return summary


def save_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2), encoding="utf-8")


def main() -> None:
    args = parse_args()
    args.experiments_root = args.experiments_root.resolve()
    evaluation_root = (
        args.evaluation_root.resolve()
        if args.evaluation_root is not None
        else default_evaluation_root(args.experiments_root)
    )

    if args.reset and evaluation_root.exists():
        shutil.rmtree(evaluation_root)
    evaluation_root.mkdir(parents=True, exist_ok=True)
    mirror_runtime_root(args.experiments_root, evaluation_root)

    lt_cmd = resolve_lt_cmd(evaluation_root, args.use_json)
    sites = available_sites(args)
    if args.max_runs is not None:
        sites = sites[: args.max_runs]
    if not sites:
        raise RuntimeError(f"No runnable site directories found under {args.experiments_root}")

    if args.watch == "interval":
        best_watch_key = "best_interval_cc"
    elif args.watch == "complement":
        best_watch_key = "best_complement_cc"
    else:
        best_watch_key = f"best_{args.watch}_cc"

    print(
        f"Running expect-style LTE supervisor over {len(sites)} site(s); "
        f"watch={args.watch} threshold={args.threshold:.3f} action={args.action} "
        f"mirror={evaluation_root}",
        flush=True,
    )

    history = []
    best_run = None
    best_value = float("-inf")

    for index, site in enumerate(sites, start=1):
        site_id = str(site["id"])
        site_name = str(site.get("name", site_id))
        print(
            f"[{index}/{len(sites)}] {site_id} - {site_name}",
            flush=True,
        )
        summary = run_site(args, lt_cmd, site, args.experiments_root, evaluation_root)
        history.append(summary)

        current_best = summary.get(best_watch_key)
        if current_best is None:
            threshold_event = summary.get("threshold_event") or {}
            current_best = threshold_event.get("watched_cc")
        if current_best is not None and current_best > best_value:
            best_value = float(current_best)
            best_run = summary

        print(
            f" completed site={site_id} exit={summary['exit_code']} "
            f"best_{metric_names(args.exclude)[0]}={summary.get(f'best_{metric_names(args.exclude)[0]}_cc')} "
            f"best_{metric_names(args.exclude)[1]}={summary.get(f'best_{metric_names(args.exclude)[1]}_cc')}",
            flush=True,
        )

        save_json(evaluation_root / "super_history.json", history)
        save_json(
            evaluation_root / "super_summary.json",
            {
                "experiments_root": str(args.experiments_root),
                "evaluation_root": str(evaluation_root),
                "results_html": str(args.results_html.resolve()) if args.results_html.exists() else None,
                "lt_cmd": lt_cmd,
                "mirror_strategy": "copy top-level runtime files and each source site directory into evaluation_root before running",
                "processed_run_count": len(history),
                "watch": args.watch,
                "threshold": args.threshold,
                "action": args.action,
                "best_run": best_run,
                "history": history,
            },
        )

    if best_run is not None:
        print(
            f"Best watched run: site={best_run['site_id']} "
            f"value={best_value:.6f}",
            flush=True,
        )


if __name__ == "__main__":
    main()
