#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Standalone TG-220-style Water-Equivalent Diameter (WED) calculator with Tkinter GUI

What it does (similar workflow to the iDose WED web page):
- Select one or more DICOM CT images (single-slice files).
- Compute per-image WED from HU + pixel spacing using TG-220-style water-equivalent area:
    Aw = sum_over_mask ( (1 + HU/1000) * pixel_area )
    WED = 2 * sqrt(Aw/pi)
- Show a list of slices with WED, allow selecting a slice and exporting results to CSV.

Dependencies:
- numpy
- pydicom
Optional (faster/cleaner morphology):
- scipy (scipy.ndimage)

Install:
  pip install numpy pydicom
Optional:
  pip install scipy
"""

from __future__ import annotations

import csv
import math
import os
import threading
import traceback
from dataclasses import dataclass
from datetime import datetime
from typing import List, Optional, Tuple

import numpy as np
import pydicom

try:
    from scipy import ndimage as ndi  # optional
    HAVE_SCIPY = True
except Exception:
    HAVE_SCIPY = False

import tkinter as tk
from tkinter import filedialog, messagebox, ttk


# -----------------------------
# Core computation
# -----------------------------

@dataclass
class WEDResult:
    path: str
    wed_cm: float
    aw_cm2: float
    rows: int
    cols: int
    pixel_spacing_mm: Tuple[float, float]
    slice_location: Optional[float]
    image_position_z: Optional[float]
    instance_number: Optional[int]
    series_description: Optional[str]
    study_date: Optional[str]


def _safe_float(x) -> Optional[float]:
    try:
        if x is None:
            return None
        return float(x)
    except Exception:
        return None


def _safe_int(x) -> Optional[int]:
    try:
        if x is None:
            return None
        return int(x)
    except Exception:
        return None


def dicom_to_hu(ds: pydicom.Dataset) -> np.ndarray:
    """
    Convert DICOM pixel data to HU using RescaleSlope/RescaleIntercept.
    """
    arr = ds.pixel_array.astype(np.float32)

    slope = float(getattr(ds, "RescaleSlope", 1.0))
    intercept = float(getattr(ds, "RescaleIntercept", 0.0))

    hu = arr * slope + intercept
    return hu


def get_pixel_spacing_mm(ds: pydicom.Dataset) -> Tuple[float, float]:
    """
    Return (row_spacing_mm, col_spacing_mm). Defaults to (1,1) if missing.
    """
    ps = getattr(ds, "PixelSpacing", None)
    if ps is None or len(ps) != 2:
        return (1.0, 1.0)
    return (float(ps[0]), float(ps[1]))


def largest_connected_component(mask: np.ndarray) -> np.ndarray:
    """
    Keep only the largest connected component in a boolean mask.
    Uses scipy.ndimage if available; otherwise a pure numpy/stack BFS fallback.
    """
    mask = mask.astype(bool)

    if mask.sum() == 0:
        return mask

    if HAVE_SCIPY:
        labeled, n = ndi.label(mask)
        if n <= 1:
            return mask
        counts = np.bincount(labeled.ravel())
        counts[0] = 0  # background
        keep = counts.argmax()
        return labeled == keep

    # Fallback (no scipy): BFS over 4-neighborhood. Slower but works for typical CT sizes.
    h, w = mask.shape
    visited = np.zeros_like(mask, dtype=bool)

    def neighbors(r, c):
        if r > 0: yield (r - 1, c)
        if r < h - 1: yield (r + 1, c)
        if c > 0: yield (r, c - 1)
        if c < w - 1: yield (r, c + 1)

    best_comp = []
    for r in range(h):
        # quick skip rows with no foreground
        if not mask[r].any():
            continue
        for c in range(w):
            if not mask[r, c] or visited[r, c]:
                continue
            stack = [(r, c)]
            visited[r, c] = True
            comp = []
            while stack:
                rr, cc = stack.pop()
                comp.append((rr, cc))
                for nr, nc in neighbors(rr, cc):
                    if mask[nr, nc] and not visited[nr, nc]:
                        visited[nr, nc] = True
                        stack.append((nr, nc))
            if len(comp) > len(best_comp):
                best_comp = comp

    out = np.zeros_like(mask, dtype=bool)
    for rr, cc in best_comp:
        out[rr, cc] = True
    return out


def make_body_mask(hu: np.ndarray, threshold_hu: float = -350.0) -> np.ndarray:
    """
    Estimate patient/body contour mask.
    Typical threshold choices:
      -350 HU: often includes lung partially; more permissive than -500
      -500 HU: stricter; may exclude very low-density regions

    Steps:
    - threshold
    - fill holes (if scipy)
    - keep largest connected component
    - optional small closing to smooth (if scipy)
    """
    # basic threshold
    mask = hu > float(threshold_hu)

    # Keep largest object to reject external artifacts/labels
    mask = largest_connected_component(mask)

    if HAVE_SCIPY:
        # Fill holes inside patient region (helps for lung cavities etc.)
        mask = ndi.binary_fill_holes(mask)

        # light morphological closing to smooth edges (structure ~ 3x3)
        structure = np.ones((3, 3), dtype=bool)
        mask = ndi.binary_closing(mask, structure=structure, iterations=1)

        # again keep largest in case closing created small bits
        mask = largest_connected_component(mask)

    return mask.astype(bool)


def compute_wed_from_dicom(path: str, threshold_hu: float = -350.0) -> WEDResult:
    """
    Compute WED for one DICOM CT image file.
    """
    ds = pydicom.dcmread(path, force=True)

    # Pixel data check
    if "PixelData" not in ds:
        raise ValueError("No PixelData in file (not an image DICOM?).")

    hu = dicom_to_hu(ds)
    rps_mm, cps_mm = get_pixel_spacing_mm(ds)
    pixel_area_cm2 = (rps_mm * cps_mm) / 100.0  # mm^2 -> cm^2 (since 1 cm^2 = 100 mm^2)

    mask = make_body_mask(hu, threshold_hu=threshold_hu)

    # TG-220-like water-equivalent area
    # (1 + HU/1000) approximates relative attenuation vs water in diagnostic energy range.
    rel = 1.0 + (hu / 1000.0)
    rel = np.clip(rel, 0.0, None)  # avoid negative contributions from very low HU regions

    aw_cm2 = float((rel[mask].sum()) * pixel_area_cm2)
    wed_cm = float(2.0 * math.sqrt(max(aw_cm2, 0.0) / math.pi))

    slice_location = _safe_float(getattr(ds, "SliceLocation", None))
    ipp = getattr(ds, "ImagePositionPatient", None)
    image_position_z = None
    if ipp is not None and len(ipp) >= 3:
        image_position_z = _safe_float(ipp[2])

    instance_number = _safe_int(getattr(ds, "InstanceNumber", None))
    series_description = getattr(ds, "SeriesDescription", None)
    study_date = getattr(ds, "StudyDate", None)

    return WEDResult(
        path=path,
        wed_cm=wed_cm,
        aw_cm2=aw_cm2,
        rows=int(hu.shape[0]),
        cols=int(hu.shape[1]),
        pixel_spacing_mm=(rps_mm, cps_mm),
        slice_location=slice_location,
        image_position_z=image_position_z,
        instance_number=instance_number,
        series_description=str(series_description) if series_description is not None else None,
        study_date=str(study_date) if study_date is not None else None,
    )


# -----------------------------
# Tkinter GUI
# -----------------------------

class WEDApp(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("TG-220 WED Calculator - Ver 1.0 - John Stratakis")
        self.geometry("1100x650")
        self.minsize(1000, 600)

        self.results: List[WEDResult] = []
        self._worker: Optional[threading.Thread] = None
        self._stop_flag = threading.Event()

        self._build_ui()

    def _build_ui(self):
        pad = {"padx": 10, "pady": 6}

        top = ttk.Frame(self)
        top.pack(fill="x", **pad)

        ttk.Label(
            top,
            text="Automated Calculation of Water-equivalent Diameter (WED) based on TG-220 (John Stratakis)",
            font=("Segoe UI", 12, "bold"),
        ).pack(anchor="w")

        # Controls
        controls = ttk.Frame(self)
        controls.pack(fill="x", **pad)

        self.threshold_var = tk.DoubleVar(value=-350.0)

        ttk.Label(controls, text="Body threshold (HU):").grid(row=0, column=0, sticky="w")
        thr = ttk.Entry(controls, textvariable=self.threshold_var, width=10)
        thr.grid(row=0, column=1, sticky="w", padx=(6, 16))

        ttk.Button(controls, text="Select DICOM files…", command=self.on_select_files).grid(
            row=0, column=2, sticky="w", padx=(0, 10)
        )
        ttk.Button(controls, text="Clear", command=self.on_clear).grid(row=0, column=3, sticky="w", padx=(0, 10))
        ttk.Button(controls, text="Export CSV…", command=self.on_export_csv).grid(row=0, column=4, sticky="w")

        controls.columnconfigure(5, weight=1)

        self.status_var = tk.StringVar(value="Ready.")
        ttk.Label(controls, textvariable=self.status_var).grid(row=0, column=6, sticky="e")

        # Progress
        prog = ttk.Frame(self)
        prog.pack(fill="x", **pad)
        self.progress = ttk.Progressbar(prog, mode="determinate")
        self.progress.pack(fill="x")

        # Table
        mid = ttk.Frame(self)
        mid.pack(fill="both", expand=True, **pad)

        columns = (
            "wed_cm",
            "aw_cm2",
            "instance",
            "slice_loc",
            "ipp_z",
            "rows",
            "cols",
            "ps_row",
            "ps_col",
            "series",
            "path",
        )

        self.tree = ttk.Treeview(mid, columns=columns, show="headings", height=18)
        self.tree.pack(side="left", fill="both", expand=True)

        headings = {
            "wed_cm": "WED (cm)",
            "aw_cm2": "Aw (cm²)",
            "instance": "Instance",
            "slice_loc": "SliceLocation",
            "ipp_z": "IPP z",
            "rows": "Rows",
            "cols": "Cols",
            "ps_row": "PixelSp Row (mm)",
            "ps_col": "PixelSp Col (mm)",
            "series": "SeriesDescription",
            "path": "File",
        }
        widths = {
            "wed_cm": 90,
            "aw_cm2": 90,
            "instance": 70,
            "slice_loc": 90,
            "ipp_z": 90,
            "rows": 60,
            "cols": 60,
            "ps_row": 120,
            "ps_col": 120,
            "series": 220,
            "path": 420,
        }

        for c in columns:
            self.tree.heading(c, text=headings[c])
            self.tree.column(c, width=widths[c], anchor="w")

        vsb = ttk.Scrollbar(mid, orient="vertical", command=self.tree.yview)
        vsb.pack(side="right", fill="y")
        self.tree.configure(yscrollcommand=vsb.set)

        # Bottom actions
        bottom = ttk.Frame(self)
        bottom.pack(fill="x", **pad)

        ttk.Button(bottom, text="Copy selected WED", command=self.on_copy_selected).pack(side="left")
        ttk.Button(bottom, text="Copy all WEDs", command=self.on_copy_all).pack(side="left", padx=(10, 0))

        self.selected_var = tk.StringVar(value="Selected: (none)")
        ttk.Label(bottom, textvariable=self.selected_var).pack(side="right")

        self.tree.bind("<<TreeviewSelect>>", self.on_select_row)

        # Info
        info = ttk.Label(
            self,
            text=(
                "Tip: If your mask excludes lung/low-density regions too aggressively, lower the threshold (e.g., -500). "
                "If it includes too much background/noise, increase it (e.g., -200)."
            ),
            foreground="#444",
        )
        info.pack(anchor="w", padx=12, pady=(0, 8))

    def on_select_files(self):
        paths = filedialog.askopenfilenames(
            title="Select DICOM CT images (single-slice files)",
            filetypes=[
                ("DICOM files", "*.dcm *.dicom *.*"),
                ("All files", "*.*"),
            ],
        )
        if not paths:
            return

        if self._worker and self._worker.is_alive():
            messagebox.showwarning("Busy", "A computation is already running.")
            return

        self._stop_flag.clear()
        self.results = []
        for item in self.tree.get_children():
            self.tree.delete(item)

        threshold = float(self.threshold_var.get())
        self.progress["value"] = 0
        self.progress["maximum"] = len(paths)
        self.status_var.set(f"Computing WED for {len(paths)} file(s)…")

        def worker():
            errors = 0
            computed: List[WEDResult] = []
            for i, p in enumerate(paths, start=1):
                if self._stop_flag.is_set():
                    break
                try:
                    r = compute_wed_from_dicom(p, threshold_hu=threshold)
                    computed.append(r)
                except Exception:
                    errors += 1
                    # keep going, but store nothing
                self._ui_progress(i)

            # Sort in a reasonable way: InstanceNumber then SliceLocation then IPP z then filename
            def sort_key(r: WEDResult):
                return (
                    r.instance_number if r.instance_number is not None else 10**9,
                    r.slice_location if r.slice_location is not None else 10**9,
                    r.image_position_z if r.image_position_z is not None else 10**9,
                    os.path.basename(r.path).lower(),
                )

            computed.sort(key=sort_key)
            self.results = computed
            self._ui_done(errors)

        self._worker = threading.Thread(target=worker, daemon=True)
        self._worker.start()

    def _ui_progress(self, i: int):
        def _():
            self.progress["value"] = i
        self.after(0, _)

    def _ui_done(self, errors: int):
        def _():
            for r in self.results:
                self.tree.insert(
                    "",
                    "end",
                    values=(
                        f"{r.wed_cm:.2f}",
                        f"{r.aw_cm2:.1f}",
                        "" if r.instance_number is None else r.instance_number,
                        "" if r.slice_location is None else f"{r.slice_location:.2f}",
                        "" if r.image_position_z is None else f"{r.image_position_z:.2f}",
                        r.rows,
                        r.cols,
                        f"{r.pixel_spacing_mm[0]:.3f}",
                        f"{r.pixel_spacing_mm[1]:.3f}",
                        "" if r.series_description is None else r.series_description,
                        r.path,
                    ),
                )

            msg = f"Done. {len(self.results)} slice(s) computed."
            if errors:
                msg += f" {errors} file(s) failed (possibly non-image or compressed DICOM)."
            self.status_var.set(msg)

        self.after(0, _)

    def on_clear(self):
        self._stop_flag.set()
        self.results = []
        for item in self.tree.get_children():
            self.tree.delete(item)
        self.progress["value"] = 0
        self.status_var.set("Cleared. Ready.")
        self.selected_var.set("Selected: (none)")

    def on_export_csv(self):
        if not self.results:
            messagebox.showinfo("No data", "No results to export.")
            return

        default_name = f"wed_results_{datetime.now().strftime('%Y%m%d_%H%M%S')}.csv"
        out = filedialog.asksaveasfilename(
            title="Save CSV",
            defaultextension=".csv",
            initialfile=default_name,
            filetypes=[("CSV", "*.csv")],
        )
        if not out:
            return

        try:
            with open(out, "w", newline="", encoding="utf-8") as f:
                w = csv.writer(f)
                w.writerow([
                    "file",
                    "wed_cm",
                    "aw_cm2",
                    "instance_number",
                    "slice_location",
                    "image_position_z",
                    "rows",
                    "cols",
                    "pixel_spacing_row_mm",
                    "pixel_spacing_col_mm",
                    "series_description",
                    "study_date",
                ])
                for r in self.results:
                    w.writerow([
                        r.path,
                        f"{r.wed_cm:.4f}",
                        f"{r.aw_cm2:.4f}",
                        "" if r.instance_number is None else r.instance_number,
                        "" if r.slice_location is None else f"{r.slice_location:.6f}",
                        "" if r.image_position_z is None else f"{r.image_position_z:.6f}",
                        r.rows,
                        r.cols,
                        f"{r.pixel_spacing_mm[0]:.6f}",
                        f"{r.pixel_spacing_mm[1]:.6f}",
                        "" if r.series_description is None else r.series_description,
                        "" if r.study_date is None else r.study_date,
                    ])
            messagebox.showinfo("Saved", f"CSV exported:\n{out}")
        except Exception as e:
            messagebox.showerror("Export failed", f"Could not write CSV:\n{e}")

    def on_select_row(self, _event=None):
        sel = self.tree.selection()
        if not sel:
            self.selected_var.set("Selected: (none)")
            return
        item = self.tree.item(sel[0])
        vals = item.get("values", [])
        if len(vals) >= 1:
            self.selected_var.set(f"Selected WED: {vals[0]} cm")

    def on_copy_selected(self):
        sel = self.tree.selection()
        if not sel:
            messagebox.showinfo("Nothing selected", "Select a row first.")
            return
        vals = self.tree.item(sel[0]).get("values", [])
        if not vals:
            return
        text = f"WED (cm): {vals[0]}\nFile: {vals[-1]}"
        self.clipboard_clear()
        self.clipboard_append(text)
        self.status_var.set("Copied selected WED to clipboard.")

    def on_copy_all(self):
        if not self.results:
            messagebox.showinfo("No data", "No results to copy.")
            return
        lines = ["wed_cm,file"]
        for r in self.results:
            lines.append(f"{r.wed_cm:.4f},{r.path}")
        self.clipboard_clear()
        self.clipboard_append("\n".join(lines))
        self.status_var.set("Copied all WEDs to clipboard.")


def main():
    try:
        app = WEDApp()
        app.mainloop()
    except Exception:
        traceback.print_exc()
        messagebox.showerror("Fatal error", traceback.format_exc())


if __name__ == "__main__":
    main()
