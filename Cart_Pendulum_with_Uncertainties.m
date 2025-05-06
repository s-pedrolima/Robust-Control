%% Cart-Pendulum with Unceartanties

close all;
clearvars;
clc;
import Controllers.*

%% System matrices
F = [1, -4.9e-6, 2.5e-3, 0;
     -1.2e-6, 1, 1.6e-4, 0;
      9.9e-5, 0, 1, 0;
      0, -9.9e-5, 0, 1];

G = [-7.82e-4;
      3.53e-4;
      0;
      0];

%% Uncertainty matrices
H_k = [0.209; 0.104; 0.018; 0.009];

E_F = [-0.1217, -0.031, -0.17, 0.0087];

E_G = 0.0029;

%% Covariance matrices
[n_size, m_size] = size(G);

P = 1e3 * eye(n_size);
Q = eye(n_size);
R = eye(m_size);

R_Bryson = 160;
Q_Bryson = 1e5 * [1 0 0 0;
                  0 0.01 0 0;
                  0 0 1 0;
                  0 0 0 0.01];

P_Bryson = 1e8 * Q_Bryson;

%% Simulation setup
Ts = 0.0001;
N = 60 / Ts; % Number of steps

% State vector
X_unc_gc = zeros(n_size, N+1);
X_unc_rlqr_Bryson = zeros(n_size, N+1);

% Control input
u_unc_rlqr = zeros(m_size, N);
u_unc_rlqr_Bryson = zeros(m_size, N);

% Position reference
Amplitude = 2 * (0.35);
Signal_period = 30;

ref_signal = Amplitude * (-gensig('square', Signal_period, N*Ts, Ts) + 0.5);

%% Figure configuration
% Figure position and size
Figure_config = struct('PaperPositionMode', 'auto', ...
                       'Units', 'centimeters', ...
                       'Position', [8 4 40 22]);

% Font size and line width
LabelFontSizeValue = 20;
AxesFontSizeValue = 16;
Object_config = struct('LineWidth', 2.5);

%% Main

% Uncertainty calculation
Delta_i = 1;
DeltaF = H_k * Delta_i * E_F;
DeltaG = H_k * Delta_i * E_G;

F = 1*(F - eye(n_size)) + eye(n_size); % Slight adjustment

% Uncertainties in F and G
F_unc = F + DeltaF * Ts;
G_unc = G + DeltaG * Ts;


% Backward calculation
epsilon = 0.000012; 
for i = 1:1000
    % RLQR 
    [L_Bryson, K_Bryson, P_Bryson] = robust_lqr_reduced(F, G, E_F, E_G, Q_Bryson, R_Bryson, P_Bryson);

    % Guaranteed Cost 
    [K, P] = guaranteed_cost(F, G, E_F, E_G, H_k, Q, R, P, epsilon);
end

% Main loop
for i = 1:N
    Position_ref = ref_signal(i) * [0; 0; 0; Ts];
    
    u_unc_rlqr(:, i) = K * X_unc_gc(:, i);
    u_unc_rlqr_Bryson(:, i) = K_Bryson * X_unc_rlqr_Bryson(:, i);
    
    X_unc_gc(:, i+1) = F_unc * X_unc_gc(:, i) + G_unc * u_unc_rlqr(:, i) + Position_ref;
    X_unc_rlqr_Bryson(:, i+1) = F_unc * X_unc_rlqr_Bryson(:, i) + G_unc * u_unc_rlqr_Bryson(:, i) + Position_ref;
end

%% Results
% Angular position
states_fig = figure();
set(states_fig, Figure_config);

subplot(3,1,1)
x1_obj = plot((1:N+1)*Ts, X_unc_gc(1,:), ...
              (1:N+1)*Ts, X_unc_rlqr_Bryson(1,:), '-.');
set(x1_obj, Object_config);
ylabel('$rad$', 'Interpreter', 'Latex', 'fontsize', LabelFontSizeValue)
xlabel('t (s)', 'Interpreter', 'Latex', 'fontsize', LabelFontSizeValue)
axis tight
ax = gca;
grid on
ax.FontSize = AxesFontSizeValue;
title('Angular position')

% Linear position
subplot(3,1,2)
x2_obj = plot((1:N+1)*Ts, X_unc_gc(2,:), ...
              (1:N+1)*Ts, X_unc_rlqr_Bryson(2,:), '-.', ...
              (1:N+1)*Ts, ref_signal, 'k:');
set(x2_obj, Object_config);
ylabel('$m$', 'Interpreter', 'Latex', 'fontsize', LabelFontSizeValue)
xlabel('t (s)', 'Interpreter', 'Latex', 'fontsize', LabelFontSizeValue)
axis tight
ax = gca;
grid on
ax.FontSize = AxesFontSizeValue;
legend('Guaranteed Cost', 'RLQR Bryson''s rule', 'Position ref', 'AutoUpdate', 'off')
title('Linear position')

% Control input
subplot(3,1,3)
u1_plot_obj = plot((1:N)*Ts, u_unc_rlqr(1,:), ...
                   (1:N)*Ts, u_unc_rlqr_Bryson(1,:), '-.');
set(u1_plot_obj, Object_config);
title('Control input')
ylabel('$u$', 'Interpreter', 'Latex', 'fontsize', LabelFontSizeValue)
xlabel('t (s)', 'Interpreter', 'Latex', 'fontsize', LabelFontSizeValue)
axis tight
ax = gca;
grid on
ax.FontSize = AxesFontSizeValue;

% Uncomment if you want to save the figure
% print('Results_RLQR_Pendulum','-dpng')


