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


## Performance Considerations

GPU execution significantly outperforms serial CPU execution for large-scale 3D simulations. However, MATLAB performance optimization remains essential. To maximize speed:

- Avoid unnecessary `for` loops  
- Minimize array indexing inside update equations  
- Prefer vectorized operations  
- Preallocate arrays before simulation  

Proper coding practices combined with GPU acceleration result in substantial computational speedup.

---

## Requirements

- MATLAB (R20XX or later recommended)
- Parallel Computing Toolbox
- CUDA-enabled NVIDIA GPU

---

---

## 📌 Citation

If you use this code in your research, please cite the following paper:

**Plain Text Citation**

Author(s),  
"Application of GPU-Accelerated FDTD Method to Electromagnetic Wave Propagation in Plasma Using MATLAB Parallel Processing Toolbox,"  
*Journal Name*, vol. XX, no. XX, pp. XX–XX, Year.  
DOI: xx.xxxx/xxxxx

---

**BibTeX**

```bibtex
@article{yourkey2024,
  author  = {Author1 and Author2 and Author3},
  title   = {Application of GPU-Accelerated FDTD Method to Electromagnetic Wave Propagation in Plasma Using MATLAB Parallel Processing Toolbox},
  journal = {Journal Name},
  volume  = {XX},
  number  = {XX},
  pages   = {XX--XX},
  year    = {2024},
  doi     = {xx.xxxx/xxxxx}
}

