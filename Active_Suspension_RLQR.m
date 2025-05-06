%% Active Suspension with Unceartanties

close all; 
clearvars; 
clc;
import Controllers.*

%% System parameters
M_s = 2.45;   % Sprung Mass
M_us = 1;     % Unsprung Mass
K_s = 900;    % Suspension Stiffness
K_us = 1250;  % Tire Stiffness
B_s = 7.5;    % Suspension Damping Coefficient
B_us = 5;     % Tire Damping Coefficient
Ts = 0.01;   % Sampling time

%% System matrices
A = [0   1     0    -1;
    -K_s/M_s  -B_s/M_s   0  B_s/M_s;
     0   0     0     1;
     K_s/M_us B_s/M_us  -K_us/M_us -(B_s+B_us)/M_us];

B = [0 0;
     0 1/M_s;
    -1 0;
     B_us/M_us -1/M_us];

[n_size, m_size] = size(B);

% Discretization
F = eye(n_size) + A * Ts;
G = B * Ts;

%% Parametric uncertainties

% Uncertainty matrices
H = [1; 1; 1; 1];
E_F = [-3.5 0.75 0.5 -0.001];
E_G = 0.01;

Deltai = 0.25*(2 * rand() - 1); 
F_unc = F+H*Deltai*E_F;
G_unc = G+H*Deltai*E_G;



%% Weighting matrices
P_r = eye(n_size);
P_g = eye(n_size);

Q = [450  0    0    0;
       0  30   0    0;
       0   0   5    0;
       0   0   0  0.01];

R = 0.01;

%% Simulation setup
N = 10 / Ts;  % Number of steps

% State vectors
X_unc_gc = zeros(n_size, N+1);
X_unc_rlqr  = zeros(n_size, N+1);

% Control inputs
u_unc_gc = zeros(m_size-1, N);
u_unc_rlqr  = zeros(m_size-1, N);

% Position reference
Amplitude = 2 * (0.01);
Signal_period = 3;
ref_signal = Amplitude * (-gensig('square', Signal_period, N*Ts, Ts) + 1);
ref_signal = [zeros(200,1); ref_signal(1:end-200)];
Road_profile = diff(ref_signal') / Ts;  % Reference is given in velocity

%% Figure configuration
Figure_config = struct('PaperPositionMode', 'auto', ...
                       'Units', 'centimeters', ...
                       'Position', [3 3 20 15]);

LabelFontSizeValue = 25;
AxesFontSizeValue  = 20;

Object_config = struct('LineWidth', 2.5);

%% Main

% Backward calculation
epsilon = 0.000000127;
for i = 1:1000
    [L_r, K_r, P_r] = robust_lqr(F, G(:,2), E_F, E_G, H, Q, R, P_r); % RLQR
    [K_g, P_g] = guaranteed_cost(F, G(:,2), E_F, E_G, H, Q, R, P_g, epsilon); % Guaranteed Cost
end

for i = 1:N
    u_unc_gc(:,i) = K_g * X_unc_gc(:,i);
    u_unc_rlqr(:,i)  = K_r * X_unc_rlqr(:,i);
    
    X_unc_gc(:,i+1) = F_unc * X_unc_gc(:,i) + G_unc * [Road_profile(i); u_unc_gc(:,i)];
    X_unc_rlqr(:,i+1)  = F_unc * X_unc_rlqr(:,i)  + G_unc * [Road_profile(i); u_unc_rlqr(:,i)];
end

%% Results

% States comparison - Top plate position
x1_fig = figure(); 
set(x1_fig, Figure_config);
x1_obj = plot((1:N+1)*Ts, X_unc_rlqr(1,:) + ref_signal' + X_unc_rlqr(3,:), '-', ...
              (1:N+1)*Ts, X_unc_gc(1,:) + ref_signal' + X_unc_gc(3,:), '--',...
              (1:N+1)*Ts, ref_signal, 'k:');
set(x1_obj, Object_config);
ylabel('$z_s(m)$', 'Interpreter', 'Latex', 'fontsize', LabelFontSizeValue)
xlabel('Time (s)', 'Interpreter', 'Latex', 'fontsize', LabelFontSizeValue)
legend('Uncertain (RLQR)', 'Uncertain (GC)')
axis tight
ax = gca;
ax.YGrid = 'on';
ax.FontSize = AxesFontSizeValue;
title('Top plate position')
% print('top_plate_pos','-dpng')

% Vehicle body acceleration
x2_diff_fig = figure();
set(x2_diff_fig, Figure_config);
x2_diff_obj = plot((1:N)*Ts, diff(X_unc_rlqr(2,:))/Ts, '-', ...
                   (1:N)*Ts, diff(X_unc_gc(2,:))/Ts, '--');
set(x2_diff_obj, Object_config);
ylabel('$m/s^2$', 'Interpreter', 'Latex', 'fontsize', LabelFontSizeValue)
xlabel('Time (s)', 'Interpreter', 'Latex', 'fontsize', LabelFontSizeValue)
legend('Uncertain (RLQR)', 'Uncertain (GC)')
axis tight
ax = gca;
ax.YGrid = 'on';
ax.FontSize = AxesFontSizeValue;
title('Vehicle body acceleration')
% print('vehicle_accel','-dpng')

% Control inputs - Active suspension force
Control_fig = figure();
set(Control_fig, Figure_config);
u1_plot_obj = plot((1:N)*Ts, u_unc_rlqr(1,:), '-', ...
                   (1:N)*Ts, u_unc_gc(1,:), '--');
set(u1_plot_obj, Object_config);
title('Active suspension control')
ylabel('$F_c(N)$', 'Interpreter', 'Latex', 'fontsize', LabelFontSizeValue)
xlabel('Time (s)', 'Interpreter', 'Latex', 'fontsize', LabelFontSizeValue)
legend('Uncertain (RLQR)', 'Uncertain (GC)')
axis tight
ax = gca;
ax.YGrid = 'on';
ax.FontSize = AxesFontSizeValue;
% print('active_suspension_control','-dpng')

% rlqr_values = eig(F_unc + G_unc(:,2)*K_r);
% gc_values = eig(F_unc + G_unc(:,2)*K_g);
% 
% Eigenvalues_fig = figure;
% set(Eigenvalues_fig, Figure_config);
% th = 0:pi/100:2*pi;
% plot(cos(th), sin(th), 'LineWidth', 2);
% hold on;
% plot(real(rlqr_values), imag(rlqr_values), '.', ...
%     'MarkerSize', 25, 'LineWidth', 2);
% plot(real(gc_values), imag(gc_values), '.', ...
%     'MarkerSize', 12, 'LineWidth', 2);
% grid on;
% axis equal;
% legend('', 'Open-loop', 'Closed-loop','Location', 'Southeast')
% title('Eigenvalues')