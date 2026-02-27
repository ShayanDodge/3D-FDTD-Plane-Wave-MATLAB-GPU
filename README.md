# 3D GPU-Accelerated FDTD in MATLAB

This repository provides a parallel, GPU-accelerated implementation of the three-dimensional Finite-Difference Time-Domain (FDTD) method using MATLAB and the Parallel Computing Toolbox.

The FDTD method is a widely used numerical technique for solving Maxwell’s equations in the time domain and is extensively applied in computational electromagnetics for simulating electromagnetic wave propagation.

---

## Overview

This implementation leverages GPU computing to significantly improve simulation performance compared to conventional CPU-based implementations. The code is optimized for high-performance execution in MATLAB while maintaining clarity and flexibility for research and educational use.

---

## Features

- **GPU Acceleration**  
  Utilizes MATLAB's Parallel Computing Toolbox to execute FDTD field updates on the GPU, significantly reducing simulation time.

- **Optimized MATLAB Implementation**  
  Code is structured to improve performance by:
  - Minimizing the use of `for` loops  
  - Reducing costly array indexing operations  
  - Exploiting vectorized computations  

- **Full 3D FDTD Solver**  
  Implements a complete three-dimensional Yee-grid formulation for electromagnetic wave propagation.

- **Modular and Readable Code**  
  Well-documented scripts designed for easy modification and research extension.

- **Visualization Tools**  
  Includes field visualization routines for monitoring electric and magnetic field evolution.

---

## Requirements

- MATLAB (R2020 or later recommended)
- Parallel Computing Toolbox
- CUDA-enabled NVIDIA GPU

---


## 📌 Citation

If you use this code in your research, please cite the following paper:

* Dodge S, Shafiee M, Shokri B. "Application of GPU-Accelerated FDTD Method to Electromagnetic Wave Propagation in Plasma Using MATLAB Parallel Processing Toolbox," arXiv preprint arXiv:2211.05647. 2022 Nov 10. DOI: [10.48550/arXiv.2211.05647](https://doi.org/10.48550/arXiv.2211.05647)

---


