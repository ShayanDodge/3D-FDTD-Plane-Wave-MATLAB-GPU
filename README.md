# FDTD 3D GPU

* This repository contains an accelerated parallel implementation of the 3D Finite Difference Time Domain (FDTD) method using GPU computing in MATLAB. The FDTD method is a powerful numerical technique for solving Maxwell's equations in the time domain, widely used in computational electromagnetics for simulating wave propagation.

## Features
* GPU Acceleration: Leverages the computational power of GPUs to significantly speed up 3D FDTD simulations.
* Optimized MATLAB Code: Includes best practices for MATLAB GPU programming, such as minimizing the use of for-loops and avoiding array indexing, to enhance execution speed.
* 3D FDTD Simulation: Implements a full three-dimensional FDTD algorithm to simulate electromagnetic wave propagation.
* User-Friendly MATLAB Scripts: Provides clear and well-documented MATLAB code, making it easy to understand and modify.
* Visualization Tools: Offers functions for visualizing electric and magnetic fields, enabling detailed analysis of wave behavior over time.
* Customizable Parameters: Allows users to adjust various simulation parameters such as grid size, time steps, and source properties.

* Although the execution speed of serial code in the CPU is slower than that of parallel code in the GPU, pay attention to some points in MATLAB programming such as avoiding For loop and avoiding array indexing improves the code execution speed.
* 
![GPU_2a](https://user-images.githubusercontent.com/94797491/152185619-5c4a3f6c-2ddb-45cf-b055-41bbe892e0e4.jpg)
