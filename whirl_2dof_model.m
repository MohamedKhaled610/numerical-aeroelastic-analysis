clear; clc; close all;


% Datum Parameter Values
R = 0.152; % Rotor Radius in m
omega_Mars = 40; % Angular Velocity in rad/s (if on Mars then = 8*omega on Earth)
a = 0.25; % Nacelle Length to Rotor Radius Ratio
I_x = 0.000103; % Rotor Moment of Inertia in kg*m^2
I_n = 0.000178; % Nacelle Moment of Inertia in kg*m^2
C_theta_fixed = 0.001; % Structural Pitch Damping in N*m*s/rad
K_theta_fixed = 0.4; % Structural Pitch Stiffness in N*m/rad
C_psi_fixed = 0.001; % Structural Yaw Damping in N*m*s/rad
K_psi_fixed = 0.4; % Structural Pitch Stiffness in N*m/rad
N_blades = 4; % Number of Blades
c = 0.026; % Blade Chord Length in m
cl_alpha = 2*pi; % Blade Lift Coefficient in 1/rad
rho = 1.225; % Density Value in kg/m^3


% Assign Velocity Range
V_vec = 0:0.2:72; 
num_steps = length(V_vec);


% Preallocate memory to store the real parts of the 4 eigenvalues
Real_Eig = zeros(4, num_steps);
Imag_Eig = zeros(4, num_steps);


% Left Side of the Equation of Motion
M_Matrix = [I_n 0; 0 I_n]; % Mass Matrix

C_Matrix = [C_theta_fixed -I_x*omega_Mars; I_x*omega_Mars C_psi_fixed]; % Damping Matrix

K_Matrix = [K_theta_fixed 0; 0 K_psi_fixed]; % Stiffness Matrix

% Derived Equations
K_a = 0.5*rho*cl_alpha*(R^4)*(omega_Mars^2);
Q = 0.5 * N_blades * K_a * R; % Common multiplier (This outputs a 2x1 array: [Min; Max])

% For State Space Matrix
Zero_Matrix = zeros(2);
I_Matrix = eye(2);

for k = 1:num_steps
    
    % Update the advance ratio for this step
    V_current = V_vec(k);
    mu = V_current / (omega_Mars*R);
    
    
    % Aerodynamic Coefficients
    integrand_A1 = @(eta) (mu^2) ./ sqrt(mu^2 + eta.^2);
    A1 = (c/R) * integral(integrand_A1, 0, 1); 
    A1_prime = mu * A1;
    
    integrand_A2_prime = @(eta) (mu^2)*(eta.^2) ./ sqrt(mu^2 + eta.^2);
    A2_prime = (c/R) * integral(integrand_A2_prime, 0, 1); 
    
    integrand_A3 = @(eta) (eta.^4) ./ sqrt(mu^2 + eta.^2);
    A3 = (c/R) * integral(integrand_A3, 0, 1); 
    
    % Aerodynamic Moments Decomposition
        
    C_diag   = Q * (-(A3 + (a^2)*A1) / omega_Mars); 
    K_direct = Q * (a * A1_prime); 
    K_cross  = Q * A2_prime; 
    
    
    % The 2x2 Matrices for Martian Density
    C_Aero = [C_diag 0; 0 C_diag];
    
    K_Aero = [K_direct -K_cross; K_cross K_direct];
    
        
    Total_C = C_Matrix - C_Aero; % Taking the Aerodynamic values to the left side of the EOM 
    Total_K = K_Matrix - K_Aero;
    
    % State Space Matrix for Martian Density
    A = [Zero_Matrix I_Matrix; -inv(M_Matrix)*Total_K -inv(M_Matrix)*Total_C];
    
    
    % Extract and Store the Imaginary and Real Parts of the Eigenvalues
    lambda = eig(A);
    Real_Eig(:, k) = real(lambda);
    Imag_Eig(:, k) = imag(lambda);

end

omega_n = sqrt(Real_Eig.^2 + Imag_Eig.^2); % Calculating the Natural Frequencies and Damping Ratios
omega_n_Hz = omega_n ./ (2*pi);
Zeta = -Real_Eig ./ omega_n;


%% Stability Boundary Analysis
V_Base = 6.7;
C_multiplier_Base = 1;
K_multiplier_Base = 1;
% Function defined at the bottom
% These are for Stiffness Diagrams
K_Grid_Baseline      = get_stiffness_envelope(V_Base,  N_blades,  C_multiplier_Base,   R,   c,  a ,     omega_Mars, rho, I_n, I_x, cl_alpha);
K_Grid_HigherDamping = get_stiffness_envelope(V_Base,  N_blades,         1.2       ,   R,   c,  a ,     omega_Mars, rho, I_n, I_x, cl_alpha);
K_Grid_HighSpeed     = get_stiffness_envelope(1.2*V_Base    ,  N_blades,  C_multiplier_Base,   R,   c,  a ,     omega_Mars, rho, I_n, I_x, cl_alpha);
K_Grid_MoreBlades    = get_stiffness_envelope(V_Base,     6    ,  C_multiplier_Base,   R,   c,  a ,     omega_Mars, rho, I_n, I_x, cl_alpha);
K_Grid_HigherOmega   = get_stiffness_envelope(V_Base,  N_blades,  C_multiplier_Base,   R,   c,  a , 1.2*omega_Mars, rho, I_n, I_x, cl_alpha);
K_Grid_Higher_a      = get_stiffness_envelope(V_Base,  N_blades,  C_multiplier_Base,   R,   c, 3*a,     omega_Mars, rho, I_n, I_x, cl_alpha);

% These are for Damping 
C_Grid_Baseline      = get_damping_envelope(V_Base,  N_blades,  K_multiplier_Base,   R,   c,  a ,     omega_Mars, rho, I_n, I_x, cl_alpha);
C_Grid_HigherStiffness = get_damping_envelope(V_Base,  N_blades,         1.2       ,   R,   c,  a ,     omega_Mars, rho, I_n, I_x, cl_alpha);
C_Grid_HighSpeed     = get_damping_envelope(1.2*V_Base   ,  N_blades,  K_multiplier_Base,   R,   c,  a ,     omega_Mars, rho, I_n, I_x, cl_alpha);C_Grid_MoreBlades    = get_damping_envelope(V_Base,     6    ,  K_multiplier_Base,   R,   c,  a ,     omega_Mars, rho, I_n, I_x, cl_alpha);
C_Grid_HigherOmega   = get_damping_envelope(V_Base,  N_blades,  K_multiplier_Base,   R,   c,  a , 1.2*omega_Mars, rho, I_n, I_x, cl_alpha);
C_Grid_Higher_a      = get_damping_envelope(V_Base,  N_blades,  K_multiplier_Base,   R,   c, 3*a,     omega_Mars, rho, I_n, I_x, cl_alpha);

%% Time Response Function
% Calculate the time response for the system using the state-space representation
V_TRF = 7.764; 

mu_TRF = V_TRF / (omega_Mars*R);


% Aerodynamic Coefficients  
integrand_A1_TRF = @(eta) (mu_TRF^2) ./ sqrt(mu_TRF^2 + eta.^2);
A1_TRF = (c/R) * integral(integrand_A1_TRF, 0, 1); 
A1_prime_TRF = mu_TRF * A1_TRF;

integrand_A2_prime_TRF = @(eta) (mu_TRF^2)*(eta.^2) ./ sqrt(mu_TRF^2 + eta.^2);
A2_prime_TRF = (c/R) * integral(integrand_A2_prime_TRF, 0, 1); 

integrand_A3_TRF = @(eta) (eta.^4) ./ sqrt(mu_TRF^2 + eta.^2);
A3_TRF = (c/R) * integral(integrand_A3_TRF, 0, 1); 

% Aerodynamic Moments Decomposition

C_diag_TRF   = Q * (-(A3_TRF + (a^2)*A1_TRF) / omega_Mars); 
K_direct_TRF = Q * (a * A1_prime_TRF); 
K_cross_TRF  = Q * A2_prime_TRF; 


% The 2x2 Matrices for Martian Density
C_Aero_TRF = [C_diag_TRF 0; 0 C_diag_TRF];

K_Aero_TRF = [K_direct_TRF -K_cross_TRF; K_cross_TRF K_direct_TRF];


Total_C_TRF = C_Matrix - C_Aero_TRF; % Taking the Aerodynamic values to the left side of the EOM 
Total_K_TRF = K_Matrix - K_Aero_TRF;

% State Space Matrix for Martian Density
A_TRF = [Zero_Matrix I_Matrix; -inv(M_Matrix)*Total_K_TRF -inv(M_Matrix)*Total_C_TRF];


t_span = 0:0.01:10; % Time vector
initial_conditions = [0.5; 0; 0; 0]; % Initial state [theta; psi; theta_dot; psi_dot]

[time_response, state_response] = ode45(@(t, y) A_TRF * y, t_span, initial_conditions);




%% Plotting

% Figure 1: Stability Boundary (Real Part of Eigenvalue)
figure('Name', 'Whirl Flutter Stability Boundary');
p_plot = plot(V_vec, Real_Eig, 'LineWidth', 2);
yline(0, 'r--', 'Flutter Boundary', 'LineWidth', 2, 'LabelHorizontalAlignment', 'left');
grid on;
xlabel('Freestream Velocity, V (m/s)');
ylabel('Damping (Real Part of Eigenvalue, \sigma)');
title('Stability at Martian Density (\rho = 0.020)');
ylim([-8 3]);
legend([p_plot(2), p_plot(4)], 'Forward Whirl Mode', 'Backward Whirl Mode', 'Location', 'southwest');

% --- Figure 2: Modal Parameters (Frequency and Damping) ---
figure('Name', 'Modal Parameters');

% Calculate the Non-Dimensional X-Axis (V / V_tip)
V_ratio = V_vec ./ (omega_Mars * R);

% Extract exact colors from Figure 1 to keep them consistent
color_fw = p_plot(2).Color;
color_bw = p_plot(4).Color;

% Subplot 1: Natural Frequency (Dual Axis) 
subplot(2, 2, 2);

% Left Y-Axis (Forward Whirl)
yyaxis left;
p_freq_fw = plot(V_ratio, omega_n_Hz(2, :), 'LineWidth', 2, 'Color', color_fw);
ylabel('Forward Whirl [Hz]');
% Set tight limits around 30 Hz to show the curvature
ylim([min(omega_n_Hz(2,:))-1, max(omega_n_Hz(2,:))+1]); 
set(gca, 'ycolor', color_fw); % Colors the axis text to match the line

% Right Y-Axis (Backward Whirl)
yyaxis right;
p_freq_bw = plot(V_ratio, omega_n_Hz(4, :), 'LineWidth', 2, 'Color', color_bw);
ylabel('Backward Whirl [Hz]');
% Set tight limits around 2 Hz to show the curvature
ylim([min(omega_n_Hz(4,:))-0.5, max(omega_n_Hz(4,:))+0.5]); 
set(gca, 'ycolor', color_bw); % Colors the axis text to match the line

grid on;
xlabel('Advance Ratio (V / V_{tip})');
title('Natural Frequency vs. Advance Ratio');


% Subplot 2: Modal Damping Ratio
subplot(2, 2, 1);
p_damp_fw = plot(V_ratio, Zeta(2, :), 'LineWidth', 2, 'Color', color_fw); hold on;
p_damp_bw = plot(V_ratio, Zeta(4, :), 'LineWidth', 2, 'Color', color_bw);

% Add a zero-line because Zeta crossing 0 is the exact flutter boundary
yline(0, 'r--', 'Flutter Boundary', 'LineWidth', 2, 'LabelHorizontalAlignment', 'left');

grid on;
xlabel('Advance Ratio (V / V_{tip})');
ylabel('Modal Damping Ratio [\zeta]');
title('Modal Damping Ratio vs. Advance Ratio');
legend([p_damp_fw, p_damp_bw], 'Forward Whirl', 'Backward Whirl', 'Location', 'southwest');
% Optional: lock the Y-axis to match your reference image scaling
ylim([-0.1 0.1]);


% Subplot 3: Stability Boundary
subplot(2, 2, 3);

stiffness_range = 0:0.005:0.5;

% Plot the 2D grid
% X-axis maps to columns (K_psi), Y-axis maps to rows (K_theta)
imagesc(stiffness_range, stiffness_range, K_Grid_Baseline);

% IMPORTANT: imagesc draws matrices from top-to-bottom by default. 
% This command flips the Y-axis so 0 is at the bottom left corner.
set(gca, 'YDir', 'normal'); 

% Set a custom two-color map: Dark Blue for 0, Bright Yellow for 1
colormap([0.1 0.2 0.5; 1 0.9 0]); 

% Setup the colorbar to act as a legend
c = colorbar;
c.Ticks = [0.25 0.75];
c.TickLabels = {'Flutter (Unstable)', 'Safe (Stable)'};

% Add Labels and Title
xlabel('Yaw Stiffness, K_\psi (N*m/rad)');
ylabel('Pitch Stiffness, K_\theta (N*m/rad)');
title(sprintf('Whirl Flutter Stability Region at V = %.1f m/s', V_Base));




% Subplot 4: Time Response Function
subplot(2,2,4);
% Y_out(:,1) is Pitch (theta), Y_out(:,2) is Yaw (psi)
plot(time_response, state_response(:, 1), 'LineWidth', 1.5); hold on; % Light Blue Solid
plot(time_response, state_response(:, 2), 'LineStyle', '--', 'LineWidth', 1.5); % Orange Dashed

grid on;
xlabel('Time [sec]', 'FontSize', 11);
ylabel('Displacement [rad]', 'FontSize', 11);
title(sprintf('Time Response \n V = %.1f m/sec', V_TRF), 'FontWeight', 'bold');
legend('\theta', '\psi', 'Location', 'northeast');
xlim([0 10]);



% --- Figure 3 PARAMETRIC COMPARISON SCRIPT (STIFFNESS ENVELOPE) ---

figure('Name', 'Parametric Stability Comparison -- Stiffness Diagrams');
hold on;
grid on;

% Plot each boundary using contour
% The [0.5 0.5] tells MATLAB to draw the line exactly between 0 and 1
contour(stiffness_range, stiffness_range, K_Grid_Baseline,       [0.5 0.5], 'y',   'LineWidth', 4); % Thick solid black
contour(stiffness_range, stiffness_range, K_Grid_HighSpeed,      [0.5 0.5], 'b--', 'LineWidth', 2.5); % Blue dashed
contour(stiffness_range, stiffness_range, K_Grid_MoreBlades,     [0.5 0.5], 'r-.', 'LineWidth', 2.5); % Red dash-dot
contour(stiffness_range, stiffness_range, K_Grid_HigherDamping,  [0.5 0.5], 'g:',  'LineWidth', 2.5); % Green dotted
contour(stiffness_range, stiffness_range, K_Grid_Higher_a,       [0.5 0.5], 'w-' ,'LineWidth', 2.5); % Green dotted
contour(stiffness_range, stiffness_range, K_Grid_HigherOmega ,   [0.5 0.5], 'o:' ,'LineWidth', 2.5); % Green dotted


% Format the axes
xlabel('Yaw Stiffness, K_\psi (N*m/rad)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Pitch Stiffness, K_\theta (N*m/rad)', 'FontSize', 12, 'FontWeight', 'bold');
title('Whirl Flutter Stability Boundaries (Terrestrial Density)', 'FontSize', 14);

% Add the Legend
legend(sprintf('Baseline (V=%.1f, 4 Blades)',V_Base), ...
    'V=1.2*V_{Base}', ...
    '6 Blades', ...
    '1.2* Structural Damping', ...
    '3*a', ...
    '1.2*omega',...
    'Location', 'best', 'FontSize', 11);

% Ensure 0 is at the bottom left
set(gca, 'YDir', 'normal'); 
hold off;



% --- Figure 4 PARAMETRIC COMPARISON SCRIPT (Damping ENVELOPE) ---

figure('Name', 'Parametric Stability Comparison -- Damping Diagrams');
hold on;
grid on;

damping_range = 0:0.00002:0.002;

% Plot each boundary using contour
% The [0.5 0.5] tells MATLAB to draw the line exactly between 0 and 1
contour(damping_range, damping_range, C_Grid_Baseline,       [0.5 0.5], 'y--',   'LineWidth', 4); % Thick solid yellow
contour(damping_range, damping_range, C_Grid_HighSpeed,      [0.5 0.5], 'b--', 'LineWidth', 2.5); % Blue dashed
contour(damping_range, damping_range, C_Grid_MoreBlades,     [0.5 0.5], 'r-.', 'LineWidth', 2.5); % Red dash-dot
contour(damping_range, damping_range, C_Grid_HigherStiffness,  [0.5 0.5], 'g:',  'LineWidth', 2.5); % Green dotted
contour(damping_range, damping_range, C_Grid_Higher_a,       [0.5 0.5], 'w' ,'LineWidth', 2.5);  
contour(damping_range, damping_range, C_Grid_HigherOmega ,   [0.5 0.5], 'o:' ,'LineWidth', 2.5); 


% Format the axes
xlabel('Yaw Damping, C_\psi (N*m s/rad)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Pitch Damping, C_\theta (N*m s/rad)', 'FontSize', 12, 'FontWeight', 'bold');
title('Whirl Flutter Stability Boundaries (Terrestrial Density)', 'FontSize', 14);

% Add the Legend
legend(sprintf('Baseline (V=%.1f, 4 Blades)',V_Base), ...
    'V=1.2*V_{Base}', ...
    '6 Blades', ...
    '1.2* Structural Stiffness', ...
    '3*a', ...
    '1.2*omega',...
    'Location', 'best', 'FontSize', 11);

% Ensure 0 is at the bottom left
set(gca, 'YDir', 'normal'); 
hold off;





figure('Name', 'Stiffness Diagram');
stiffness_range = 0:0.005:0.5;

% Plot the 2D grid
% X-axis maps to columns (K_psi), Y-axis maps to rows (K_theta)
imagesc(stiffness_range, stiffness_range, K_Grid_Baseline);

% IMPORTANT: imagesc draws matrices from top-to-bottom by default. 
% This command flips the Y-axis so 0 is at the bottom left corner.
set(gca, 'YDir', 'normal'); 

% Set a custom two-color map: Dark Blue for 0, Bright Yellow for 1
colormap([0.1 0.2 0.5; 1 0.9 0]); 

% Setup the colorbar to act as a legend
c = colorbar;
c.Ticks = [0.25 0.75];
c.TickLabels = {'Flutter (Unstable)', 'Safe (Stable)'};

% Add Labels and Title
xlabel('Yaw Stiffness, K_\psi (N*m/rad)');
ylabel('Pitch Stiffness, K_\theta (N*m/rad)');
title(sprintf('Whirl Flutter Stability Region at V = %.1f m/s', V_Base));


figure('Name', 'Time Response');
% Y_out(:,1) is Pitch (theta), Y_out(:,2) is Yaw (psi)
plot(time_response, state_response(:, 1), 'LineWidth', 1.5); hold on; % Light Blue Solid
plot(time_response, state_response(:, 2), 'LineStyle', '--', 'LineWidth', 1.5); % Orange Dashed

grid on;
xlabel('Time [sec]', 'FontSize', 11);
ylabel('Displacement [rad]', 'FontSize', 11);
title(sprintf('Time Response \n V = %.3f m/sec', V_TRF), 'FontWeight', 'bold');
legend('\theta', '\psi', 'Location', 'northeast');
xlim([0 10]);



%% Exact Flutter Velocity Extraction (Physical Modes)
fprintf('\n--- Whirl Flutter Stability Report ---\n');

% Execute the extraction for both Martian density cases
[V_critical, V_all] = extract_critical_V(Real_Eig, V_vec);

% Print Results for Density (rho = 0.020)
fprintf('\n--- rho = %.3f ---\n', rho);
print_physical_modes(V_all, V_critical);



% --- Helper Function 1: Console Formatting ---
function print_physical_modes(V_all, V_critical)
    % Group the conjugate pairs (1 & 2 -> Physical Mode 1, 3 & 4 -> Physical Mode 2)
    physical_mode_1 = V_all(1); 
    physical_mode_2 = V_all(3); 
    
    if ~isnan(physical_mode_1)
        fprintf('Physical Mode 1 (Forward Whirl) : Flutter at %.3f m/s\n', physical_mode_1);
    else
        fprintf('Physical Mode 1 (Forward Whirl)  : Stable across all tested speeds.\n');
    end
    
    if ~isnan(physical_mode_2)
        fprintf('Physical Mode 2 (Backward Whirl) : Flutter at %.3f m/s\n', physical_mode_2);
    else
        fprintf('Physical Mode 2 (Backward Whirl)  : Stable across all tested speeds.\n');
    end
    
    if ~isnan(V_critical)
        fprintf('>> ABSOLUTE CRITICAL VELOCITY    : %.3f m/s\n', V_critical);
    end
end


% --- Helper Function 2: Zero-Crossing Extraction ---
function [V_crit_lowest, V_crit_all] = extract_critical_V(Eig_Matrix, V_vec)
    V_crit_all = NaN(1, 4); 
    
    for mode = 1:4
        % Find the first index where the real part crosses zero
        cross_idx = find(Eig_Matrix(mode, 1:end-1) < 0 & Eig_Matrix(mode, 2:end) >= 0, 1);
        
        if ~isempty(cross_idx)
            % Extract the coordinates surrounding the crossing
            v1 = V_vec(cross_idx);
            v2 = V_vec(cross_idx + 1);
            r1 = Eig_Matrix(mode, cross_idx);
            r2 = Eig_Matrix(mode, cross_idx + 1);
            
            % Linear interpolation formula for exact root
            V_exact = v1 - r1 * (v2 - v1) / (r2 - r1);
            
            % Save the exact velocity to its specific mode slot in the array
            V_crit_all(mode) = V_exact;
        end
    end
    
    % The absolute critical velocity is just the lowest number in the array
    V_crit_lowest = min(V_crit_all);
end




function [Results_Diagram] = get_stiffness_envelope(V_static, num_blades, C_multiplier, R_rotor, c, a, omega, rho, I_n, I_x, cl)

C_theta_fixed = 0.001;
C_psi_fixed = 0.001;
Zero_Matrix = zeros(2);
I_Matrix = eye(2);

mu_static = V_static / (omega*R_rotor);
K_a = 0.5*rho*cl*(R_rotor^4)*(omega^2);
Q = 0.5 * num_blades * K_a * R_rotor;


% Aerodynamic Coefficients
integrand_A1_static = @(eta) (mu_static^2) ./ sqrt(mu_static^2 + eta.^2);
A1_static = (c/R_rotor) * integral(integrand_A1_static, 0, 1); 
A1_prime_static = mu_static * A1_static;

integrand_A2_prime_static = @(eta) (mu_static^2)*(eta.^2) ./ sqrt(mu_static^2 + eta.^2);
A2_prime_static = (c/R_rotor) * integral(integrand_A2_prime_static, 0, 1); 

integrand_A3_static = @(eta) (eta.^4) ./ sqrt(mu_static^2 + eta.^2);
A3_static = (c/R_rotor) * integral(integrand_A3_static, 0, 1); 

% Aerodynamic Moments Decomposition

C_diag_static   = Q * (-(A3_static + (a^2)*A1_static) / omega); 
K_direct_static = Q * (a * A1_prime_static); 
K_cross_static  = Q * A2_prime_static; 


% The 2x2 Matrices for Martian Density
C_Aero_static = [C_diag_static 0; 0 C_diag_static];

K_Aero_static = [K_direct_static -K_cross_static; K_cross_static K_direct_static];

% Left Side of the Equation of Motion
M_Matrix = [I_n 0; 0 I_n]; % Mass Matrix


C_theta = C_theta_fixed * C_multiplier;
C_psi = C_psi_fixed * C_multiplier;
C_Matrix = [C_theta -I_x*omega; I_x*omega C_psi]; % Damping Matrix


Total_C_static = C_Matrix - C_Aero_static;


num_steps= length(0:0.005:0.5);
Results_Diagram = zeros(num_steps,num_steps);

row = 1;

for K_theta = 0:0.005:0.5

    column = 1;

    for K_psi = 0:0.005:0.5


        K_Matrix_variable = [K_theta 0; 0 K_psi]; % Varying Stiffness Matrix

        Total_K_variable = K_Matrix_variable - K_Aero_static;

        % State Space Matrix for Martian Density
        A_varying_K = [Zero_Matrix I_Matrix; -inv(M_Matrix)*Total_K_variable -inv(M_Matrix)*Total_C_static];

        % Determine stability by checking whether the eigenvalues became
        % non negative or no

        Real_Eig_variable_max = max(real(eig(A_varying_K)));

        if Real_Eig_variable_max < 0
            Results_Diagram(row,column) = 1;
        else
            Results_Diagram(row,column) = 0;
        end


        column = column + 1; % Move to the next column for K_psi

    end
    row = row + 1; % Move to the next row for K_theta

end

end



function [Results_Diagram] = get_damping_envelope(V_static, num_blades, K_multiplier, R_rotor, c, a, omega, rho, I_n, I_x, cl)

K_theta_fixed = 0.4;
K_psi_fixed = 0.4;

Zero_Matrix = zeros(2);
I_Matrix = eye(2);

mu_static = V_static / (omega*R_rotor);
K_a = 0.5*rho*cl*(R_rotor^4)*(omega^2);
Q = 0.5 * num_blades * K_a * R_rotor;

num_steps= length(0:0.005:0.5);


% Aerodynamic Coefficients
integrand_A1 = @(eta) (mu_static^2) ./ sqrt(mu_static^2 + eta.^2);
A1_static = (c/R_rotor) * integral(integrand_A1, 0, 1); 
A1_prime_static = mu_static * A1_static;

integrand_A2_prime = @(eta) (mu_static^2)*(eta.^2) ./ sqrt(mu_static^2 + eta.^2);
A2_prime_static = (c/R_rotor) * integral(integrand_A2_prime, 0, 1); 

integrand_A3 = @(eta) (eta.^4) ./ sqrt(mu_static^2 + eta.^2);
A3_static = (c/R_rotor) * integral(integrand_A3, 0, 1); 

% Aerodynamic Matrices
C_diag_static   = Q * (-(A3_static + (a^2)*A1_static) / omega); 
K_direct_static = Q * (a * A1_prime_static); 
K_cross_static  = Q * A2_prime_static; 

C_Aero_static = [C_diag_static 0; 0 C_diag_static];
K_Aero_static = [K_direct_static -K_cross_static; K_cross_static K_direct_static];

% 4. Lock the Fixed Stiffness Matrix (Calculated OUTSIDE the loop)
M_Matrix = [I_n 0; 0 I_n]; 

K_theta = K_multiplier * K_theta_fixed;
K_psi = K_multiplier * K_psi_fixed;
K_Matrix = [K_theta 0; 0 K_psi];

Total_K_static = K_Matrix - K_Aero_static;

% Pre-allocate the grid
Results_Diagram = zeros(num_steps, num_steps);

row = 1;

for C_theta = 0:0.00002:0.002

    column = 1;

    for C_psi = 0:0.00002:0.002

        C_Matrix_variable = [C_theta -I_x*omega; I_x*omega C_psi]; 

        Total_C_variable = C_Matrix_variable - C_Aero_static;

        % State Space Matrix
        A_varying_C = [Zero_Matrix I_Matrix; -inv(M_Matrix)*Total_K_static -inv(M_Matrix)*Total_C_variable];

        % Determine stability
        Real_Eig_variable_max = max(real(eig(A_varying_C)));

        if Real_Eig_variable_max < 0
            Results_Diagram(row, column) = 1; % Stable
        else
            Results_Diagram(row, column) = 0; % Unstable
        end

        column = column + 1; 
    end
    row = row + 1; 
end
end