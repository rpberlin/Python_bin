#!/usr/bin/env python3
import csv
import glob
import math
import os
import shlex
import signal
import subprocess
import sys
import time
from typing import Dict, List, Optional, Tuple

import matplotlib.pyplot as plt

HISTORY_FILE = "history.csv"
XCOL = "Inner_Iter"
MONITOR_COLS = ["CD", "CL", "CSF", "CMx", "CMy", "CMz", "CFx", "CFy", "CFz", "CEff"]

def which(cmd: str) -> Optional[str]:
    for p in os.environ.get("PATH", "").split(os.pathsep):
        cand = os.path.join(p, cmd)
        if os.path.isfile(cand) and os.access(cand, os.X_OK):
            return cand
    return None

def default_su2_exe() -> str:
    env = os.environ.get("SU2_CFD_EXE")
    if env:
        return env
    if which("SU2_CFD_84"):
        return "SU2_CFD_84"
    return "SU2_CFD"

SU2_EXE = default_su2_exe()

def norm(s: str) -> str:
    return s.strip().strip('"').strip("'")

def safe_float(x: str) -> float:
    t = norm(x)
    if not t:
        return float("nan")
    t = t.replace("D", "E").replace("d", "E")
    try:
        return float(t)
    except Exception:
        return float("nan")

def is_residual_col(name: str) -> bool:
    return norm(name).lower().startswith("rms")

def pick_cfg() -> str:
    cfgs = sorted(glob.glob("*.cfg"))
    if not cfgs:
        print("No .cfg files found in current directory.")
        sys.exit(1)
    if len(cfgs) == 1:
        print(f"Found one cfg: {cfgs[0]} (auto-selected)")
        return cfgs[0]
    print("Available .cfg files:")
    for i, c in enumerate(cfgs, start=1):
        print(f"  {i:2d}) {c}")
    while True:
        s = input("Select cfg number: ").strip()
        if s.isdigit():
            idx = int(s)
            if 1 <= idx <= len(cfgs):
                return cfgs[idx - 1]
        print("Invalid selection. Try again.")

def ask_run_mode() -> Tuple[List[str], int]:
    have_mpirun = which("mpirun") is not None
    print("\nRun mode:")
    print("  1) Serial")
    print(f"  2) MPI (mpirun) {'(available)' if have_mpirun else '(NOT found on PATH)'}")
    while True:
        mode = input("Choose 1 or 2: ").strip()
        if mode == "1":
            return [], 1
        if mode == "2" and have_mpirun:
            np_ = int(input("MPI ranks (-np): ").strip())
            t_s = input("OpenMP threads per rank (-t) [default 1]: ").strip() or "1"
            t_ = int(t_s)
            return ["mpirun", "-np", str(np_)], t_
        print("Invalid selection. Try again.")

def finite(seq: List[float]) -> List[float]:
    return [v for v in seq if isinstance(v, (int, float)) and math.isfinite(v)]

class HistoryTail:
    """
    Robust tailer: reopens if file is truncated/rotated.
    """
    def __init__(self, path: str):
        self.path = path
        self.fp = None
        self.header: Optional[List[str]] = None
        self.pos = 0
        self.last_size = -1

    def _open(self) -> bool:
        if not os.path.exists(self.path):
            return False
        self.fp = open(self.path, "r", newline="")
        self.fp.seek(0)
        # read first non-empty line as header
        while True:
            line = self.fp.readline()
            if not line:
                return False
            if line.strip():
                raw = next(csv.reader([line]))
                self.header = [norm(h) for h in raw]
                self.pos = self.fp.tell()
                self.last_size = os.path.getsize(self.path)
                return True

    def open_if_ready(self) -> bool:
        if self.fp is not None and self.header is not None:
            # detect truncation/rotation
            try:
                sz = os.path.getsize(self.path)
                if sz < self.last_size:
                    self.close()
                    return self._open()
                self.last_size = sz
                return True
            except Exception:
                return True
        return self._open()

    def read_new(self) -> List[List[str]]:
        if self.fp is None or self.header is None:
            return []
        self.fp.seek(self.pos)
        rows = []
        reader = csv.reader(self.fp)
        for r in reader:
            if r and any(cell.strip() for cell in r):
                rows.append([norm(c) for c in r])
        self.pos = self.fp.tell()
        return rows

    def close(self):
        if self.fp:
            self.fp.close()
        self.fp = None
        self.header = None
        self.pos = 0
        self.last_size = -1

def main():
    cfg = pick_cfg()
    prefix, threads = ask_run_mode()

    cmd: List[str] = []
    cmd.extend(prefix)
    cmd.append(SU2_EXE)
    if prefix:
        cmd.extend(["-t", str(threads)])
    cmd.append(cfg)

    stamp = time.strftime("%Y%m%d_%H%M%S")
    log_name = f"run_{os.path.splitext(cfg)[0]}_{stamp}.log"
    print(f"Command: {cmd}")

    print("\nLaunching:")
    print("  " + " ".join(shlex.quote(x) for x in cmd))
    print(f"Logging to: {log_name}")
    print(f"Watching: {HISTORY_FILE} (x={XCOL}, residuals=rms* are LOG10 already)\n")

    logf = open(log_name, "w")
    proc = subprocess.Popen(cmd, stdout=logf, stderr=subprocess.STDOUT, text=True)


    default_w, default_h = plt.rcParams["figure.figsize"]
    fig, (ax_res, ax_mon) = plt.subplots(
        1, 2,
        figsize=(default_w * 1.8, default_h * 1.2),
        sharex=True
    )
    plt.ion()
    fig.tight_layout(pad=2.0)

    tail = HistoryTail(HISTORY_FILE)

    x: List[float] = []
    data: Dict[str, List[float]] = {}
    residual_cols: List[str] = []
    monitor_cols: List[str] = []

    last_plot = 0.0
    plot_period_s = 1.0
    max_pts = 30000

    def update_plots():
        if not x:
            return
        xs = x[-max_pts:]

        # Residuals (already log10 values, so linear axis)
        ax_res.cla()
        ax_res.set_title("Residuals")
        ax_res.set_ylabel("log10(rms[...])")
        ax_res.set_xlabel(XCOL)
        ax_res.grid(True, alpha=0.25)

        ymins, ymaxs = [], []
        for c in residual_cols:
            ys_all = data.get(c, [])
            ys = ys_all[-max_pts:] if ys_all else []
            f = finite(ys)
            if not f:
                continue
            ax_res.plot(xs, ys, label=c)
            ymins.append(min(f))
            ymaxs.append(max(f))

        if ymins:
            ax_res.legend(loc="best", fontsize="small")
            ymin, ymax = min(ymins), max(ymaxs)
            pad = 0.1 * (ymax - ymin) if ymax > ymin else 0.5
            ax_res.set_ylim(ymin - pad, ymax + pad)

        # Monitors
        ax_mon.cla()
        ax_mon.set_title("Force / Moment monitors")
        ax_mon.set_xlabel(XCOL)
        ax_mon.set_ylabel("Coefficient")
        ax_mon.grid(True, alpha=0.25)

        any_mon = False
        for c in monitor_cols:
            ys_all = data.get(c, [])
            ys = ys_all[-max_pts:] if ys_all else []
            if finite(ys):
                ax_mon.plot(xs, ys, label=c)
                any_mon = True
        if any_mon:
            ax_mon.legend(loc="best", fontsize="small")
            ax_mon.relim()
            ax_mon.autoscale_view()

        fig.canvas.draw_idle()
        plt.pause(0.75)

    try:
        while True:
            if tail.open_if_ready() and tail.header:
                if not residual_cols and not monitor_cols:
                    header = tail.header
                    hset = set(header)
                    residual_cols = [h for h in header if is_residual_col(h)]
                    monitor_cols = [c for c in MONITOR_COLS if c in hset]

                    if XCOL not in hset:
                        print(f"WARNING: '{XCOL}' not found; using row index for x-axis.")

                    for h in residual_cols + monitor_cols:
                        data.setdefault(h, [])

                    print(f"Detected residual cols: {residual_cols}")
                    print(f"Detected monitor cols:  {monitor_cols}")

                new_rows = tail.read_new()
                if new_rows and tail.header:
                    idx = {name: i for i, name in enumerate(tail.header)}
                    for r in new_rows:
                        if len(r) < len(tail.header):
                            continue
                        if XCOL in idx:
                            x.append(safe_float(r[idx[XCOL]]))
                        else:
                            x.append(float(len(x)))
                        for c in residual_cols + monitor_cols:
                            if c in idx:
                                data[c].append(safe_float(r[idx[c]]))

            now = time.time()
            if now - last_plot >= plot_period_s:
                last_plot = now
                update_plots()

            if proc.poll() is not None:
                break

            time.sleep(0.05)

    except KeyboardInterrupt:
        print("\nCtrl-C: sending SIGINT to SU2 for graceful stop...")
        try:
            proc.send_signal(signal.SIGINT)
        except Exception:
            proc.terminate()
        proc.wait()

    finally:
        tail.close()
        logf.close()
        print(f"\nRun finished (return code {proc.returncode}). Log: {log_name}")
        plt.ioff()
        plt.show()

if __name__ == "__main__":
    main()
