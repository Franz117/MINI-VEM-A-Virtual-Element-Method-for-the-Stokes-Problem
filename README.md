# MINI-VEM: A Virtual Element Method for the Stokes Problem

⚠️ **Non-commercial use only**  
This code is provided for academic and research purposes. Commercial use is strictly prohibited.

## Overview
MINI-VEM is a MATLAB implementation of a **novel Virtual Element Method** inspired by the classical Mini element for the **incompressible Stokes equations**. The method works on general polygonal meshes and ensures stable velocity-pressure coupling.

This repository contains scripts for:
- Mesh generation (triangular, quadrilateral, polygonal)
- MINI-VEM, MINI, quad MINI (Lamanchè), UNVIRTUALIZED MINI-VEM.
- Post-processing and visualization
- Validation and convergence studies

## Folder Structure
- `main/` — Main scripts to run simulations
- `src/UTIL/` — Core routines and helper functions
- `meshes/` — Mesh generation scripts and saved meshes
- `results/` — Numerical results (.mat files)
- `images/` — Figures and diagrams
- `obsolete/` — Deprecated code

## Usage
1. Open MATLAB and add the repo folder to the path:
   ```matlab
   addpath(genpath('MINI-VEM'))