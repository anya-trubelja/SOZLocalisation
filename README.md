# SOZLocalisation
MATLAB toolbox for stereo‑EEG seizure‑onset zone localisation

# StereoEEG‑Based SOZ Localization & Multimodal Neuroimaging Toolbox

A collection of MATLAB scripts and helper functions for end‑to‑end seizure‑onset‑zone (SOZ) localization using stereo‑EEG, Freesurfer surfaces, SPM, 3D Slicer, and geodesic analyses.

---

## 🔍 Overview

For a full step-by-step guide here is the user manual:
https://docs.google.com/document/d/12HZBheDz0HKaEXPKnc4RbHN7O6ibOMVZa_2d6WhbBvo/edit?usp=sharing

This toolbox supports the full pipeline for focal epilepsy surgical planning:

1. **DICOM → NIfTI Conversion**  
   Import MRI/CT into 3D Slicer and export as NIfTI.

2. **Electrode Segmentation & Export**  
   Use GARDEL to segment electrodes on CT, export MRI‑centered coordinates to Excel.

3. **Flatmap Patch Generation**  
   - Extract Desikan–Killiany (DKT) ROIs via Freesurfer.  
   - Build `*.dkt.electrode.patch` & `*.flat.patch` from selected regions.

4. **Electrode Contact Mapping**  
   - Compute nearest‐pial vertex for each contact.  
   - Export per‐contact Excel tables.

5. **Time‑Difference‑of‑Arrival (TDOA) Multilateration**  
   - Parse electrode onset times (HH:MM:SS:cc).  
   - Compute propagation distances & all 3‑contact combinations.  
   - Trilaterate 2D “origins” on the flatmap, project back to 3D.

6. **Heatmap Generation & Export**  
   - Build Gaussian‑blurred density maps on the pial mesh.  
   - Export NIfTI volumes (`heatmap_nosmooth.nii`, `heatmap_smooth.nii`) via SPM.

7. **Resection Overlap Analysis**  
   - Resample resection mask into heatmap space.  
   - Threshold heatmap (90th percentile) → predicted SOZ mask.  
   - Compute true/false overlap counts, volumes, and percentages.

8. **Quality & Comparison Tools**  
   - Compare two heatmaps (.fig) vertex‑by‑vertex (difference, Pearson R).  
   - Compute local geodesic‐to‐flatmap distortion per vertex.  
   - Compare mean absolute distortion across sessions.

---

## ⚙️ Prerequisites

- **Operating System**: macOS or Linux (Freelyfer and 3D Slicer compatibility)  
- **MATLAB** ≥ R2018b with toolboxes:
  - Image Processing
  - Signal Processing
  - Statistics & Machine Learning
  - Control System, Deep Learning, Bioinformatics (as needed)
  - EEGLAB
- **Freelyfer** (via `$FREESURFER_HOME` & `$SUBJECTS_DIR`)
- **SPM12** on MATLAB path
- **GARDEL** (Windows‑only; MATLAB R2018b + SPM8)
- **3D Slicer** (for DICOM handling & segmentation)
- **Additional MATLAB toolboxes** (wavelets, graph)

---

## 🚀 Installation

1. **Clone this repository**  
   ```bash
   git clone https://github.com/your‑org/SOZLocalizationToolbox.git
   cd SOZLocalizationToolbox

2. Set up Freesurfer
  export FREESURFER_HOME=/Applications/freesurfer
  source $FREESURFER_HOME/SetUpFreeSurfer.sh
  export SUBJECTS_DIR=$FREESURFER_HOME/subjects

3. Add MATLAB paths
  addpath(genpath('Toolboxes/spm12-master'));
  addpath(genpath('Toolboxes/freesurferMatlabLibrary-master'));
  addpath(genpath('Toolboxes/toolbox_wavelets'));
  addpath(genpath('Toolboxes/toolbox_graph'));

5. Install GARDEL (if needed) and ensure SPM8 is on your MATLAB path.

📂 Folder Structure
/SOZLocalizationToolbox
│
├─ Toolboxes/
│   ├─ spm12-master/
│   ├─ freesurferMatlabLibrary-master/
│   ├─ toolbox_wavelets/
│   └─ toolbox_graph/
│
├─ scripts/
│   ├─ convertDicomToNifti.m       # 3D Slicer‑assisted export
│   ├─ createDKTElectrodeFlatmap.m
│   ├─ computeElectrodeContactTable.m
│   ├─ ATMultilateration.m
│   ├─ makeSOZNifti.m
│   ├─ compareHeatmaps.m
│   ├─ ResectionHeatmapEvaluation.m
│   ├─ LocalGeodesicDistortion.m
│   └─ compare_mean_distortion_from_figs.m
│
├─ examples/                       # Example data & .fig files
│
└─ docs/
    └─ User Manual.docx            # Detailed step‑by‑step guide

📜 License
This project is licensed under the MIT License. See LICENSE for details.


© 2025 Anya Trubelja – Improving Neurosurgical Planning in Focal Epilepsy.



