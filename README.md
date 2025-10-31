# 🛰️ Advanced 3D Extended Kalman Filter (EKF) — UAV Tracking in XYZ Space

## 📘 Overview
This project implements an **Extended Kalman Filter (EKF)** for **3D UAV tracking** in Cartesian coordinates (X, Y, Z).  
The nonlinear UAV motion model is simulated in MATLAB, and the EKF algorithm is implemented in **C/C++ for Vivado HLS**, suitable for FPGA acceleration.

It extends previous projects:
- [Advanced-3D-Particle-Filter-UAV-Tracking-in-XYZ-Space](https://github.com/abbasfadavi/Advanced-3D-Particle-Filter-UAV-Tracking-in-XYZ-Space)
- [Advanced-3D-Kalman-Filter-UAV-Tracking-in-XYZ-Space](https://github.com/abbasfadavi/Advanced-3D-Kalman-Filter-UAV-Tracking-in-XYZ-Space)

---

## ⚙️ System Description
- **State vector (6×1):**
  \[
  x = [x, y, z, v_x, v_y, v_z]^T
  \]
- **Process model:** nonlinear UAV motion with acceleration and turning dynamics  
- **Measurement model:** direct noisy position measurement  

---


