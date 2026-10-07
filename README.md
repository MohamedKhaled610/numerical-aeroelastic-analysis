# numerical-whirl-flutter-aeroelastic-analysis
A MATLAB based aeroelastic analysis suite that evaluates linear equations of motion for single-rotor nacelles to identify whirl flutter boundaries and mode coalescence. The solver conducts parametric stiffness and damping grid sweeps using state space eigenvalue tracking, and validates using time domain response function for given initial disturbance.


## Features / Solvers
Quasi Steady Aerodynamics: Computes aerodynamic moment integrals (A_1, A'_1, A'_2, A_3) using the integral() function.
Flutter Identification: Detects flutter instability by checking the sign of the real eignevalues extracted by the eig() function.
State Space Eigenvalue Grid Sweep: Converts coupled second-order equations of motion into a first order state space form (y_dot = A * y)
Time Response Function: Employs the ode45 function that uses fourth and fifth Runge-Kutta approximation in solving the first order state space equations to simulate transient pitch-yaw response after initial disturbance.   


## Requirements
- MATLAB Version: MATLAB R2022a or newer
- Toolboxes: None (Base MATLAB only)

## How to Run
1. Clone or download the repository.
2. Open MATLAB and set the directory to this folder.
3. Alter the Datum parameters as needed or use directly.
4. Run `whirl_2dof_model.m` in the command window.
