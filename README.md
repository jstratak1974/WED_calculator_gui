# TG-220 Water-Equivalent Diameter Calculator

A standalone Python desktop application for calculating the water-equivalent
diameter, $D_w$, of individual axial CT images using an AAPM TG-220-style
method.

The program reads one or more single-frame DICOM CT files, converts the stored
pixel values to Hounsfield units, estimates the patient contour, calculates the
water-equivalent cross-sectional area, and reports the corresponding
water-equivalent diameter for each slice. Results can be reviewed in the GUI,
copied to the clipboard, or exported to CSV.

> [!CAUTION]
> This program should be validated with representative phantom and patient
> images before use in clinical QA or dose reporting. The result depends on HU
> calibration, patient segmentation, reconstruction field of view, DICOM
> metadata, and the presence of SciPy. The application calculates $D_w$ only;
> it does not calculate SSDE or patient absorbed dose.

## Contents

- [Features](#features)
- [Calculation method](#calculation-method)
- [Processing workflow](#processing-workflow)
- [Requirements](#requirements)
- [Installation](#installation)
- [Running the application](#running-the-application)
- [Using the GUI](#using-the-gui)
- [Threshold selection](#threshold-selection)
- [DICOM attributes](#dicom-attributes)
- [Results and CSV export](#results-and-csv-export)
- [Worked verification example](#worked-verification-example)
- [Using the calculation function from Python](#using-the-calculation-function-from-python)
- [Validation recommendations](#validation-recommendations)
- [Current limitations](#current-limitations)
- [Troubleshooting](#troubleshooting)
- [References](#references)

## Features

- Graphical interface built with Tkinter.
- Selection of one or more single-slice DICOM CT images.
- Conversion of stored DICOM pixel values to HU using
  `RescaleSlope` and `RescaleIntercept`.
- Configurable body-contour threshold in HU.
- Rejection of disconnected external objects by retaining the largest connected
  foreground component.
- Optional hole filling and morphological smoothing with SciPy.
- Per-slice calculation of water-equivalent area and diameter.
- Background processing so the main GUI remains responsive during a batch.
- Slice table containing WED, water-equivalent area, position, dimensions,
  pixel spacing, series description, and file path.
- CSV export with four decimal places for WED and water-equivalent area.
- Clipboard export for one selected result or all calculated WED values.
- Pure-Python connected-component fallback when SciPy is unavailable.

## Calculation method

The implementation follows the axial-image concept described in AAPM Report
220. For each pixel $i$ inside the segmented patient region, the program uses
the CT number to approximate its attenuation relative to water:

$$
\frac{\mu_i}{\mu_w} \approx 1 + \frac{HU_i}{1000}.
$$

The water-equivalent cross-sectional area is calculated as

$$
A_w = \sum_{i \in ROI}
\left(1 + \frac{HU_i}{1000}\right) A_{pixel},
$$

where:

| Symbol | Meaning | Unit |
| --- | --- | --- |
| $HU_i$ | CT number of pixel $i$ | HU |
| $A_{pixel}$ | Physical area of one image pixel | cm² |
| $ROI$ | Automatically generated patient/body mask | — |
| $A_w$ | Water-equivalent cross-sectional area | cm² |

The water-equivalent diameter is the diameter of a circle with area $A_w$:

$$
D_w = 2\sqrt{\frac{A_w}{\pi}}.
$$

The script labels this result **WED (cm)** and stores it as `wed_cm`.

### Pixel area

DICOM `PixelSpacing` contains row and column spacing in millimetres. The script
converts the product to square centimetres:

$$
A_{pixel} =
\frac{\Delta r_{mm}\,\Delta c_{mm}}{100}.
$$

For example, a pixel spacing of `0.8 × 0.8 mm` gives
$A_{pixel}=0.0064\;cm^2$.

### HU conversion

The stored pixel array is converted to `float32`, then transformed using:

$$
HU = PixelValue \times RescaleSlope + RescaleIntercept.
$$

If either rescale attribute is missing, the script uses these defaults:

| Attribute | Default |
| --- | ---: |
| `RescaleSlope` | 1.0 |
| `RescaleIntercept` | 0.0 |

Relative attenuation values below zero are clipped to zero before integration:

```python
rel = 1.0 + hu / 1000.0
rel = np.clip(rel, 0.0, None)
```

Therefore, pixels below −1000 HU do not make a negative contribution to
$A_w$.

## Processing workflow

For every selected file, the application performs these steps:

1. Read the dataset using `pydicom.dcmread(path, force=True)`.
2. Confirm that the dataset contains `PixelData`.
3. Decode `pixel_array` and convert it to HU.
4. Read row and column pixel spacing.
5. Create an initial mask using `HU > threshold`.
6. Keep the largest four-connected foreground component.
7. If SciPy is installed:
   - fill internal holes;
   - apply one 3 × 3 binary-closing iteration; and
   - keep the largest component again.
8. Calculate relative attenuation and clip negative values to zero.
9. Sum the water-equivalent pixel areas inside the mask.
10. Convert $A_w$ to $D_w$.
11. Extract optional slice and series metadata.
12. Sort all successful results and display them in the table.

The batch is sorted by:

1. `InstanceNumber`;
2. `SliceLocation`;
3. the z component of `ImagePositionPatient`; and
4. filename.

Missing values are placed after available values for the corresponding sort
field.

## Requirements

- Python 3.9 or later is recommended.
- NumPy.
- pydicom.
- Tk/Tkinter.
- SciPy is optional in the code, but recommended for consistent body masking.

The basic dependencies are:

```text
numpy
pydicom
scipy
```

The application uses the Python standard library for its GUI, CSV export,
threading, and data structures.

## Installation

Clone or download the repository, then create a virtual environment:

```bash
python -m venv .venv
```

Activate the environment on Windows PowerShell:

```powershell
.venv\Scripts\Activate.ps1
```

Activate it on Linux or macOS:

```bash
source .venv/bin/activate
```

Install the recommended dependencies:

```bash
python -m pip install --upgrade pip
python -m pip install numpy pydicom scipy
```

### Tkinter on Linux

Tkinter is commonly supplied with Python on Windows and macOS. On Debian or
Ubuntu, install it separately when required:

```bash
sudo apt update
sudo apt install python3-tk
```

Confirm that Tkinter works:

```bash
python -m tkinter
```

A small Tk test window should appear.

### Compressed DICOM images

Native uncompressed transfer syntaxes can normally be decoded using pydicom and
NumPy alone. JPEG, JPEG-LS, JPEG 2000, and some other compressed transfer
syntaxes may require an additional decoder.

One common pydicom configuration is:

```bash
python -m pip install pylibjpeg pylibjpeg-libjpeg pylibjpeg-openjpeg pylibjpeg-rle
```

Install only the decoders required and approved for the local environment, and
verify decoded pixel values against a trusted DICOM viewer or reference
implementation.

## Running the application

From the directory containing the script:

```bash
python wed-new.py
```

On systems where Python 3 is invoked explicitly:

```bash
python3 wed-new.py
```

The main window opens with a default body threshold of `−350 HU`.

## Using the GUI

1. Enter the desired **Body threshold (HU)**.
2. Select **Select DICOM files…**.
3. Choose one or more axial, single-frame CT image files from the same series.
4. Wait for the progress bar and status message to show completion.
5. Review the per-slice values in the table.
6. Select a row to display its WED at the lower right.
7. Use **Copy selected WED**, **Copy all WEDs**, or **Export CSV…**.

Selecting a new group of files clears the previous batch. The current GUI does
not append new files to existing results.

### Buttons

| Button | Action |
| --- | --- |
| **Select DICOM files…** | Selects and processes a new batch of files. |
| **Clear** | Requests the worker to stop, clears current results, resets progress, and clears the selected-value label. |
| **Export CSV…** | Writes the full result table to a timestamped CSV file. |
| **Copy selected WED** | Copies the selected WED and its source file path. |
| **Copy all WEDs** | Copies two-column, CSV-like text containing WED and file path. |

## Threshold selection

The initial body mask uses:

```python
mask = hu > threshold_hu
```

The comparison is strictly greater than the threshold. A pixel exactly equal
to the threshold is excluded from the initial mask.

The default is `−350 HU`.

| Threshold change | Typical effect on initial mask |
| --- | --- |
| More negative, for example −500 HU | Includes more low-density anatomy and may include more external material or noise. |
| Less negative, for example −200 HU | Excludes more low-density anatomy and may remove valid lung or peripheral tissue. |

The final effect also depends on whether SciPy is available.

### With SciPy

The application fills holes inside the largest component and applies a small
morphological closing. Enclosed low-density regions can therefore remain within
the patient ROI and contribute according to their HU values. This is especially
relevant in thoracic images.

### Without SciPy

The application performs thresholding and largest-component selection only.
There is no hole filling or morphological closing. Low-density anatomy that
falls below the threshold may remain excluded, which can change $A_w$ and
$D_w$ relative to a SciPy-enabled installation.

For reproducible datasets, keep the dependency environment and SciPy version
consistent and record the selected threshold.

### Table and external-object removal

Keeping the largest connected component often removes detached couch
structures, labels, and external artifacts. It does not explicitly segment the
CT table. An object connected to the patient contour can remain in the mask,
and an unusually large external component could be selected instead of the
patient. Visual mask review is therefore needed during validation.

## DICOM attributes

The following attributes affect the calculation or exported metadata.

| DICOM attribute | Tag | Required by script | Use |
| --- | --- | --- | --- |
| `PixelData` | (7FE0,0010) | Yes | Source encoded pixel data. |
| `Rows` / `Columns` | (0028,0010) / (0028,0011) | Indirectly | The result records the decoded array dimensions. |
| `PixelSpacing` | (0028,0030) | No | Row and column spacing used to calculate physical pixel area. |
| `RescaleIntercept` | (0028,1052) | No | Converts stored values to HU; defaults to 0.0. |
| `RescaleSlope` | (0028,1053) | No | Converts stored values to HU; defaults to 1.0. |
| `InstanceNumber` | (0020,0013) | No | Display, CSV export, and primary sorting. |
| `SliceLocation` | (0020,1041) | No | Display, CSV export, and secondary sorting. |
| `ImagePositionPatient` | (0020,0032) | No | Its third value is displayed/exported as `image_position_z` and used for tertiary sorting. |
| `SeriesDescription` | (0008,103E) | No | Display and CSV export. |
| `StudyDate` | (0008,0020) | No | CSV export. |

### Missing pixel spacing

If `PixelSpacing` is absent or does not have exactly two values, the function
silently substitutes `(1.0, 1.0) mm`:

```python
if ps is None or len(ps) != 2:
    return (1.0, 1.0)
```

This fallback can produce a substantially incorrect physical area and WED.
Confirm `PixelSpacing` before accepting any result. A production version should
reject missing or invalid spacing rather than substituting a nominal value.

### File acceptance

The reader uses `force=True` and checks only for `PixelData`. It does not verify
the CT Image Storage SOP Class, the `Modality` value, series consistency, or
single-frame status. Select only known axial CT images from a consistent series.

## Results and CSV export

### GUI columns

| Column | Description |
| --- | --- |
| `WED (cm)` | Calculated water-equivalent diameter, displayed to two decimal places. |
| `Aw (cm²)` | Water-equivalent area, displayed to one decimal place. |
| `Instance` | DICOM `InstanceNumber`, when present. |
| `SliceLocation` | DICOM `SliceLocation`, when present. |
| `IPP z` | Third component of `ImagePositionPatient`, when present. |
| `Rows` / `Cols` | Decoded image-array dimensions. |
| `PixelSp Row (mm)` | First `PixelSpacing` value. |
| `PixelSp Col (mm)` | Second `PixelSpacing` value. |
| `SeriesDescription` | DICOM series description, when present. |
| `File` | Full source file path. |

### CSV columns

The exported CSV contains:

```text
file
wed_cm
aw_cm2
instance_number
slice_location
image_position_z
rows
cols
pixel_spacing_row_mm
pixel_spacing_col_mm
series_description
study_date
```

Numeric precision in the CSV is:

| Field | Decimal places |
| --- | ---: |
| `wed_cm` | 4 |
| `aw_cm2` | 4 |
| `slice_location` | 6 |
| `image_position_z` | 6 |
| Pixel spacing | 6 |

The default filename follows this pattern:

```text
wed_results_YYYYMMDD_HHMMSS.csv
```

The exported source path, series description, and study date may contain
sensitive or identifying information. Apply the institution's de-identification
and data-handling requirements before sharing a CSV or clipboard output.

## Worked verification example

Consider a synthetic axial image with:

- a circular water region of radius 20 pixels;
- `PixelSpacing = [1.0, 1.0] mm`;
- water pixels at `0 HU`;
- surrounding air at `−1000 HU`; and
- the circular region fully inside the image boundaries.

The discrete circle contains 1257 pixels. Each pixel has an area of
$0.01\;cm^2$, and water has a relative attenuation of 1.0:

$$
A_w = 1257 \times 0.01 = 12.57\;cm^2.
$$

Therefore:

$$
D_w = 2\sqrt{\frac{12.57}{\pi}} = 4.00058\;cm.
$$

The script returned:

```text
Aw  = 12.569999695 cm²
WED = 4.000577544 cm
```

The small difference from the analytical value is consistent with `float32`
pixel arithmetic. This test checks the implemented area and diameter equations;
it does not validate segmentation performance on clinical images.

## Using the calculation function from Python

The filename `wed-new.py` contains a hyphen, so it cannot be imported with a
normal Python `import wed-new` statement. Rename or copy it to a valid module
name such as `wed_calculator.py` before importing it:

```python
from wed_calculator import compute_wed_from_dicom

result = compute_wed_from_dicom(
    "CT000001.dcm",
    threshold_hu=-350.0,
)

print(f"WED: {result.wed_cm:.4f} cm")
print(f"Aw:  {result.aw_cm2:.4f} cm²")
print(f"Pixel spacing: {result.pixel_spacing_mm} mm")
```

`compute_wed_from_dicom()` returns a `WEDResult` dataclass with these fields:

```text
path
wed_cm
aw_cm2
rows
cols
pixel_spacing_mm
slice_location
image_position_z
instance_number
series_description
study_date
```

The computational functions can be used without launching the GUI because
`main()` runs only when the file is executed as a script.

## Validation recommendations

Before routine use:

1. Verify HU conversion for representative files, including negative stored
   values and non-default rescale slope/intercept.
2. Confirm `PixelSpacing` and the resulting pixel area against the DICOM header.
3. Test uniform water cylinders of known diameter and more complex
   tissue-equivalent phantoms.
4. Compare the generated patient mask with the CT image at the head, thorax,
   abdomen, pelvis, shoulders, and slices containing external devices.
5. Evaluate threshold sensitivity for each protocol and reconstruction type.
6. Check truncated anatomy and small-field reconstructions separately.
7. Compare SciPy-enabled and fallback results; standardize one environment for
   production calculations.
8. Verify compressed DICOM decoding against a trusted implementation.
9. Compare per-slice $D_w$ with an independent TG-220 implementation.
10. Check slice ordering using `ImagePositionPatient`, especially when instance
    numbering is absent, duplicated, or inconsistent.
11. Confirm how contrast material, metal, implants, and the patient table are
    handled in the local use case.
12. Record the code revision, dependency versions, threshold, series UID, and
    validation status with reported results.

## Current limitations

- The application calculates per-slice $D_w$ only. It does not calculate SSDE,
  CTDIvol corrections, organ dose, effective dose, or scan-average $D_w$.
- There is no image or contour display for visual verification of the mask.
- The body contour is a threshold-based largest-component estimate rather than
  a full TG-220 production segmentation workflow.
- Results can differ depending on whether SciPy is installed.
- Missing `PixelSpacing` silently becomes `1 × 1 mm`.
- Missing rescale attributes silently use slope 1 and intercept 0.
- The program does not verify `Modality = CT`, SOP Class UID, Series Instance
  UID, Frame of Reference UID, image orientation, or study/series consistency.
- Enhanced multiframe CT is not supported by the two-dimensional mask code.
- Localizer images, secondary captures, dose screens, and non-CT images are not
  explicitly rejected if they contain pixel data.
- `force=True` relaxes DICOM file-format validation.
- Pixel padding values are not removed explicitly.
- The patient table is not segmented explicitly.
- Cropped patient anatomy is not detected.
- Failed files are counted, but their filenames and exception details are not
  shown in the GUI or exported.
- Sorting gives priority to `InstanceNumber`; this may not equal anatomical
  z order in every dataset.
- Selecting files from several series can mix them into one result table.
- **Copy all WEDs** generates comma-separated text without CSV quoting, so a
  path containing a comma may not parse correctly. **Export CSV…** uses
  Python's `csv.writer` and handles quoting correctly.
- Clearing a running batch is cooperative; already queued progress or completion
  callbacks may still update the interface.
- No calculation log, mask image, dependency record, or audit file is saved.

## Troubleshooting

### `ModuleNotFoundError: No module named 'numpy'`

Install NumPy in the same Python environment used to run the script:

```bash
python -m pip install numpy
```

### `ModuleNotFoundError: No module named 'pydicom'`

Install pydicom:

```bash
python -m pip install pydicom
```

### `No module named tkinter`

Install the operating system's Tk package. On Debian/Ubuntu:

```bash
sudo apt install python3-tk
```

### Compressed DICOM files fail

Install a decoder compatible with the transfer syntax. To inspect the transfer
syntax in Python:

```python
import pydicom

ds = pydicom.dcmread("CT000001.dcm", stop_before_pixels=True)
print(ds.file_meta.TransferSyntaxUID)
print(ds.file_meta.TransferSyntaxUID.name)
```

Consult the pydicom compressed-pixel-data documentation for the corresponding
plugin.

### WED is unexpectedly low in the thorax

Check the body mask, threshold, presence of SciPy, HU calibration, pixel
spacing, and field-of-view truncation. Without SciPy, regions below the
threshold are not restored by hole filling.

### WED is unexpectedly high

Check whether the table, immobilization equipment, pads, cables, or external
objects are connected to the patient mask. Also verify that the selected file
is an axial CT image and that pixel values were converted to HU correctly.

### The status reports failed files without details

The worker intentionally catches each exception and only increments an error
counter. For debugging, call `compute_wed_from_dicom()` directly in a Python
session or modify the exception handler to record the filename and traceback.

### Slice order is incorrect

The current sorter prioritizes `InstanceNumber`. Check that the selected files
belong to one series and inspect `InstanceNumber`, `ImagePositionPatient`, and
`ImageOrientationPatient`. For robust spatial ordering, calculate position
along the slice normal from orientation and position rather than relying on the
displayed z coordinate alone.

## References

1. McCollough C, Bakalyar DM, Bostani M, et al. *Use of Water Equivalent
   Diameter for Calculating Patient Size and Size-Specific Dose Estimates
   (SSDE) in CT.* AAPM Report No. 220. 2014.
   [AAPM report page](https://www.aapm.org/pubs/reports/detail.asp?docid=146) ·
   [Report PDF](https://www.aapm.org/pubs/reports/RPT_220.pdf) ·
   [doi:10.37206/146](https://doi.org/10.37206/146)
2. American Association of Physicists in Medicine. *Size-Specific Dose
   Estimates (SSDE) in Pediatric and Adult Body CT Examinations.* AAPM Report
   No. 204. 2011.
   [Report PDF](https://www.aapm.org/pubs/reports/RPT_204.pdf)
3. pydicom project. *Handling of compressed pixel data.*
   [pydicom documentation](https://pydicom.github.io/pydicom/stable/guides/user/image_data_handlers.html)



