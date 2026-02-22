# ECG R-Wave Detection using Matched Filter

Implementation of an ECG R-peak detection algorithm based on a matched filter approach.

## Project Overview

This project implements an R-wave detection pipeline for ECG signals using:

- Bandpass filtering (5–15 Hz)
- Signal normalization
- Template-based matched filtering
- Correlation-based detection
- Adaptive thresholding
- Peak refinement and validation

The goal is to accurately detect R-peaks while minimizing false positives and false negatives.

---

## Algorithm Pipeline

1. Bandpass filtering to remove baseline wander and high-frequency noise  
2. Signal normalization to [-1,1]  
3. Template extraction from a strong R-wave  
4. Sliding correlation (matched filter)  
5. Adaptive thresholding  
6. Duplicate and distance-based filtering  
7. Final refinement using raw ECG signal  

---

## Files

- `main_R_detection.m`  
  Main execution script.

- `Rwave_detection_M.m`  
  Matched filter implementation.

- `Rwave_detection_D.m`  
  Detection and refinement logic.

---

## Technologies

- MATLAB
- Digital Signal Processing
- Matched Filtering
- Adaptive Thresholding
- ECG Signal Analysis

---

## Key Concepts

- Bandpass filter design
- Template matching
- Correlation analysis
- False positive reduction
- Distance-based peak validation

---

## Notes

Tested on provided ECG datasets.  
Designed for educational and research purposes.
