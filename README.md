# ⚙️ System Identification & Adaptive Control

A collection of three assignments implementing parameter estimation, adaptive observers, and system identification algorithms, developed for the course *Simulation and Modeling of Dynamic Systems* (2025-2026) at the Aristotle University of Thessaloniki, Department of Electrical and Computer Engineering[cite: 17, 22].

---

## 📋 Table of Contents

- [Project Structure](#project-structure)
- [Exercise 1 — Least Squares Estimation (LSE)](#exercise-1--least-squares-estimation-lse)
- [Exercise 2 — Gradient & Lyapunov Estimators](#exercise-2--gradient--lyapunov-estimators)
- [Final Project — Robust Adaptive Observers & Model Selection](#final-project--robust-adaptive-observers--model-selection)
- [Results & Key Highlights](#results--key-highlights)
- [Installation](#installation)
- [Usage](#usage)

---

## Project Structure

```text
.
├── Exercise_1/                          # Least Squares Estimation
│   ├── figures/
│   ├── Lab01_2026.pdf
│   ├── exercise_1.m                     # Forward simulation (mass-spring-damper)
│   ├── exercise_2.m                     # Inverse problem (LSE evaluation)
│   ├── least_squares.m                  # LSE algorithm implementation
│   ├── report_1.pdf
│   ├── simulate_system.m
│   └── simulation_data.mat
│
├── Exercise_2/                          # Gradient and Lyapunov Methods
│   ├── figures/
│   ├── Lab02_2026.pdf
│   ├── exercise_1.m                     # Gradient estimator for aerial vehicle
│   ├── exercise_2.m                     # Lyapunov estimator for simple pendulum
│   ├── gradient_estimator.m
│   ├── lyapunov_estimator.m
│   ├── report_2.pdf
│   ├── simulate_pendulum.m
│   └── simulate_vehicle.m
│
├── Project/                             # Robust Control & Cross-Validation
│   ├── figures/
│   ├── Project_2026.pdf
│   ├── basis_estimator.m                # Linear-in-parameters basis estimator
│   ├── lyapunov_estimator.m             # Standard adaptive observer
│   ├── part_1.m                         # Robust estimation (Box-Projection)
│   ├── part_2.m                         # Model structure selection
│   ├── projection_estimator.m           # Projection operator implementation
│   ├── report_3.pdf
│   ├── simulate_mimo.m
│   └── simulate_unknown.m
└──
```

---

## Exercise 1 — Least Squares Estimation (LSE)

Investigates the forward and inverse modeling of a classic mechanical Mass-Spring-Damper system[cite: 2, 3].
- **Forward Simulation:** Derives the state-space representation and simulates the system's dynamic response to a sinusoidal input $u(t) = 5 \sin(2.5t)$ using `ode45`[cite: 2, 3].
- **Parameter Estimation:** Implements the Least Squares Estimation (LSE) algorithm to identify the unknown mass $m$, damping $c$, and spring constant $k$[cite: 2, 3]. 
- **Robustness Analysis:** Evaluates the estimator's performance under varying sampling periods $T_s$ and differing levels of White Gaussian Noise (from 0.1% to 10%)[cite: 2, 3]. 

---

## Exercise 2 — Gradient & Lyapunov Estimators

Focuses on real-time adaptive parameter estimation using two fundamental continuous-time adaptive laws[cite: 4].
- **Gradient Method:** Applied to the linearized orientation dynamics of an aerial vehicle[cite: 4]. Estimates unknown damping $k$ and input gain $b$[cite: 4]. The system's robustness is further analyzed under the influence of an unmodeled external disturbance $d(t)$[cite: 4].
- **Lyapunov Method:** Applied to a nonlinear simple pendulum model[cite: 4, 5]. Estimates the unknown pendulum length $l$ and damping coefficient $c$[cite: 4]. The estimator's accuracy is tested against sinusoidal measurement noise $\eta(t)$[cite: 4].

---

## Final Project — Robust Adaptive Observers & Model Selection

Explores advanced adaptive control theory and practical system identification[cite: 6, 7].

### Part A: Robust Real-Time Estimation
- **Adaptive Observer:** Designs a Series-Parallel Lyapunov observer for a nonlinear MIMO system to estimate 6 unknown parameters[cite: 6, 7].
- **Box-Projection Operator:** To counter the destructive "parameter drift" caused by bounded external disturbances, a Projection Operator is implemented[cite: 7]. This constrains parameter updates within a known geometric set $\Omega$, guaranteeing Uniform Ultimate Boundedness (UUB) and ensuring absolute system stability[cite: 7].

### Part B: Model Structure Selection & Generalization
- **System Identification:** Attempts to approximate a completely unknown nonlinear scalar function $f(x,u)$ using linear combinations of basis functions[cite: 6, 7].
- **Cross-Validation:** Evaluates four candidate structures (Linear, Polynomial, Correct Structure, Overparameterized) across distinct Training, Validation, and Test datasets[cite: 7].
- **Bias-Variance Trade-off:** Demonstrates the Parsimony Principle by analyzing underfitting in simpler models and severe overfitting in overparameterized models, ultimately selecting the structure that minimizes the modeling error on unseen data[cite: 7].

---

## Results & Key Highlights

*   **Derivative Sensitivity in LSE:** Exercise 1 demonstrated that LSE is highly vulnerable to noise because approximating acceleration via velocity differentiation acts as a high-pass filter, destroying the regression matrix when noise exceeds 1%[cite: 2].
*   **Parameter Drift vs. UUB:** In the Final Project, persistent external DC disturbances caused the plain Lyapunov estimator's parameters to diverge infinitely[cite: 7]. The Projection Operator successfully clipped these values at their theoretical boundaries, sacrificing a small amount of tracking accuracy to preserve total parametric stability[cite: 7].
*   **Overfitting in System ID:** During cross-validation, the overparameterized model (9 parameters) achieved near-zero training error ($K_{train} = 0.0009$) but suffered catastrophic failure on the validation set ($K_{val} = 243.8954$), proving that extra degrees of freedom capture noise rather than underlying physics[cite: 7].
*   **Identifiability Limits:** The optimal model accurately predicted future states (Generalization index $K_{test} / K_{val} = 1.179$), but some internal estimated parameters deviated from their theoretical values due to linear dependencies between basis functions in the restricted training space[cite: 7].

---

## Installation

### Requirements

- MATLAB (Tested on R2023a+)
- Control System Toolbox (Optional, for state-space analysis)

Clone the repository to your local machine:
```bash
git clone [https://github.com/theodoratzina/System-Identification-Adaptive-Control.git](https://github.com/theodoratzina/System-Identification-Adaptive-Control.git)
```

---

## Usage

Navigate to the respective exercise directories in MATLAB and execute the main scripts:

```matlab
% Exercise 1: Least Squares Estimation
cd Exercise_1
exercise_1  % Runs the forward simulation
exercise_2  % Runs LSE evaluation and noise analysis

% Exercise 2: Gradient and Lyapunov Methods
cd Exercise_2
exercise_1  % Aerial vehicle gradient estimation
exercise_2  % Pendulum Lyapunov estimation

% Final Project: Robust Control & System Identification
cd Project
part_1      % MIMO Adaptive Observer & Projection Operator
part_2      % Model Selection & Cross-Validation
```
---
