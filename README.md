# Cartpole Control Project

## 1. Project Overview

This project focuses on the modeling, simulation, and control of a Cart-Pole system (an inverted pendulum on a cart). The objective is to stabilize the pendulum in the upright position while tracking a desired position reference for the cart, ensuring all physical constraints (like track limits and actuator saturation) are respected.

The development follows a Model-Based Design (MBD) approach, transitioning from physical equations to nonlinear simulations, continuous-time control design, and finally to discrete-time algorithms suitable for deployment on a microcontroller.

## 2. Project Structure

The project is organized into the following directories:

* **`docs/`**: Contains theoretical documentation, derivations, and initial sketches.
    * `dynamic_scheme.pdf`: Initial hand-drawn sketches of the physical system, free-body diagrams, and derivation of the nonlinear kinematic equations using Newton-Euler mechanics.
    * `dynamic_system_representation.ipynb`: A Jupyter Notebook (using Python and SymPy/NumPy) detailing the step-by-step mathematical derivation of the system's equations of motion, linearization, and state-space representation.

* **`models/`**: Contains Simulink models for system simulation.
    * `Microcontroller.slx`: A discrete-time subsystem representing the full embedded control logic (Swing-up, Feedforward, LQI, and discrete nonlinear observer).
    * `sil.slx`: Software-in-the-Loop (SIL) model to test the discrete microcontroller code against the continuous plant simulation.

* **`scripts/`**: Contains MATLAB scripts and functions for initialization, control synthesis, and simulation.
    * `cartpole_system_init.m`: Initializes all physical parameters, computes linearized matrices, synthesizes LQI and LQE (Kalman Filter) gains (both continuous and discrete).
    * `cartpole_plant.m`: Contains the nonlinear ODEs representing the physical cartpole system, including Stribeck friction and external disturbances.
    * `controller.m`: The continuous-time control law (Swing-up + LQI + Feedforward + Observer).
    * `controller_discrete_step.m`: The discrete-time equivalent of the controller, representing the exact C-code logic for the MCU.
    * `main_simulation.m`: Runs the continuous-time nonlinear simulation using `ode15s`.
    * `main_simulation_discrete.m`: Runs the discrete-time simulation, accurately modeling the Zero-Order Hold (ZOH) behavior of the MCU.
    * `system_dynamics.m`: A wrapper function combining the plant and the controller for continuous ODE simulation.
    * `track_limit.m`: An ODE event function to detect if the cart hits the physical track limits.

## 3. Mathematical Modeling & Derivation

### 3.1 Initial Analysis (`dynamic_scheme.pdf`)

The project began with a physical analysis of the system, modeled as two rigid bodies:

1. **Cart** (Mass $M$): Moves horizontally (coordinate $x$). Subject to control force $u$ and viscous/Coulomb friction.
2. **Pendulum** (Mass $m$, Length to CoM $l$): Rotates around a pivot on the cart (angle $\theta$).

Using Newton-Euler equations (balancing forces and torques), the nonlinear equations of motion were derived. These equations relate the accelerations ($\ddot{x}$, $\ddot{\theta}$) to the states ($x, \dot{x}, \theta, \dot{\theta}$), the input ($u$), and external disturbances.

### 3.2 Python Notebook Derivation (`dynamic_system_representation.ipynb`)

The hand-derived equations were formalized and verified using SymPy in a Jupyter Notebook. The process involved:

1. **Defining Symbols**: Defining all physical parameters ($M, m, l, g$, friction coefficients) and state variables.
2. **Nonlinear Equations formulation**: Constructing the inertia matrix $M(\theta)$ and the force vector $F(\theta, \dot{\theta}, \dot{x}, u)$ such that $M(\theta)\ddot{q} = F$.
3. **Equilibrium Points**: Identifying the upright position ($\theta = 0$) as the unstable equilibrium point for stabilization.
4. **Linearization**: Using Jacobian matrices (Taylor series expansion) to linearize the nonlinear equations around the upright equilibrium point.
5. **State-Space Representation**: Extracting the continuous-time state-space matrices ($A, B, C, D$) required for linear control design.

## 4. Control Architecture

The control strategy evolved in stages to handle the system's nonlinearities and physical constraints.

### 4.1 Linear Quadratic Integral (LQI) Control

Initially, an LQR controller was designed based on the linearized state-space model. However, to ensure zero steady-state error for the cart position $x$ (especially in the presence of constant disturbances like friction), an integral state was added, upgrading the controller to an LQI design.

* **State Vector**: $x_{aug} = [x_{error}, \dot{x}, \theta, \dot{\theta}, \int x_{error} dt]^T$
* **Weight Matrices (Q, R)**: Tuned using Bryson's Rule based on the maximum allowed deviations (e.g., track limits, maximum angle, maximum control effort).

### 4.2 Constant-Gain Nonlinear Observer (Hybrid Design)

To estimate the full state vector from limited sensor measurements (typically just $x$ and $\theta$), the system implements a Constant-Gain Nonlinear Observer. This architecture was deliberately chosen over a canonical Extended Kalman Filter (EKF) to respect the computational constraints of an embedded real-time target:

* **Offline Optimal Gain Synthesis**: Stationary optimal observer gains ($L$ / $L_d$) are pre-computed offline using LQE, based on the linearized system matrices around the upright equilibrium and the noise covariances ($Q_n, R_n$).
* **Nonlinear Prediction Step**: To ensure accuracy across the global domain (including far from equilibrium during the swing-up maneuver), the observer's prediction step evaluates the exact nonlinear equations of motion, incorporating the state-dependent mass matrix $M(\theta)$. 
* **Linear Correction Step ($O(1)$ complexity)**: The correction step applies the pre-computed static linear gains to the measurement innovation. This hybrid approach avoids the $O(N^3)$ computational bottleneck of calculating Jacobians and inverting covariance matrices in real-time, ensuring a deterministic $O(1)$ execution footprint strictly compatible with RTOS deadlines.

### 4.3 Energy-Based Swing-Up and Feedforward

The LQI controller is only valid near the upright position. For global control:

* **Swing-Up Controller**: An energy-shaping control law pumps energy into the pendulum to swing it up from the downward resting position ($\theta = \pi$).
* **Feedforward Friction Compensation**: To counteract the highly nonlinear Coulomb and Stribeck friction, an estimate of the normal force and resulting friction is computed and fed forward to the control signal.
* **C-infinity Blending**: A smooth, convex blending function seamlessly transitions between the swing-up controller and the stabilizing LQI controller as the pendulum approaches the upright position.

### 4.4 Discrete Implementation (MCU Target)

To prepare the algorithms for deployment on a real microcontroller, the entire control and estimation logic was discretized:

* **Plant Discretization**: The continuous plant model was discretized using Zero-Order Hold (ZOH) at the intended sampling time $T_s$.
* **Discrete LQI (DLQR)**: Optimal gains were re-calculated using `dlqr` based on the discrete matrices ($A_d, B_d$).
* **Discrete Kalman Gain Synthesis (DLQE)**: Observer gain matrices were synthesized using `dlqe`, accounting for the discrete sampling of process and measurement noise covariances.
* **Euler Integration**: Integral action and state prediction in the microcontroller code rely on explicit Forward Euler integration for computational efficiency ($O(1)$ complexity).

## 5. File Analysis & Description

### `cartpole_system_init.m`

This initialization script sets up the workspace environment.

* **Physical Parameters**: Defines masses, lengths, and highly detailed friction models (viscous, Coulomb, Stribeck effect).
* **Linearization**: Computes the analytical $A$ and $B$ matrices derived from the Python notebook.
* **LQI Synthesis**: Uses Bryson's rule to populate $Q_{aug}$ and $R_{aug}$ matrices, then computes the continuous optimal gain $K_aug$.
* **Kalman Gain Design**: Defines process ($Q_n$) and measurement ($R_n$) noise covariances based on sensor resolution, then computes the steady-state observer gain $L$ via LQE.
* **Discrete Design**: Re-computes $K_d$ and $L_d$ for the specified sampling time $T_s$ using exact ZOH discretization.

### `cartpole_plant.m`

The nonlinear simulation engine representing physics.

* Computes the exact normal force $N$ considering centripetal and tangential accelerations.
* Applies the complex Stribeck friction model.
* Solves the mass matrix system $M\ddot{q} = F$ to find accelerations for `ode15s`.

### `controller.m`

The continuous-time implementation of the control logic.

* Normalizes angles to prevent wrap-around errors.
* Computes feedforward friction compensation.
* Evaluates both the LQI stabilizing effort and the Energy-based swing-up effort.
* Blends them using a smooth exponential weight based on the angular error.
* Applies anti-windup clamping to the integral action if actuators saturate.
* Updates the nonlinear continuous observer.

### `controller_discrete_step.m`

The core MCU algorithm. This function represents exactly what runs inside the timer interrupt of the microcontroller.

* Takes discrete samples $y_{meas}$ and previous states $x_{hat}, x_i$.
* Calculates the required control effort using the discrete gains $K_d$.
* Performs Forward Euler integration for the anti-windup protected integral state.
* Runs the prediction step of the observer using the nonlinear equations.
* Runs the correction step (innovation) using the discrete pre-computed gain $L_d$ and the newly acquired measurements.

### `main_simulation.m` & `main_simulation_discrete.m`

These scripts orchestrate the simulation loops.

* The continuous version uses `ode15s` to solve the entire coupled system (plant + observer + controller) continuously.
* The discrete version simulates a real MCU: it samples the plant, runs `controller_discrete_step` once, holds the output ($ZOH$), and then integrates the plant physics forward by $T_s$ seconds. This explicitly models the sampling delays and discrete dynamics.
* Both generate comprehensive plots showing cart position (vs track limits), pendulum angle, and actuator effort (with saturation limits).

## 6. How to Run

1. Open MATLAB and ensure the `scripts/` directory is on your path.
2. Run `cartpole_system_init.m` to load all parameters into the workspace.
3. To run the continuous ideal simulation, execute `main_simulation.m`.
4. To run the high-fidelity discrete microcontroller simulation, execute `main_simulation_discrete.m`.
5. To view the mathematical derivations, open `dynamic_system_representation.ipynb` in Jupyter Notebook or VS Code.
6. To run the Simulink models, open `sil.slx` after running the initialization script.