% Robust Linear Quadratic Regulator (RLQR) Example

close all;
clearvars;
clc;
import Controllers.*

NUSP = 1+0+7+4+8+1+9+8; % ID USP

%% System Configuration 1 

% System matrices
F = [1.1, 0, 0;
     0, 0, 1.2;
     -1, 1, 0];

G = [0, 1;
     1, 1;
    -1, 1];

% Uncertainty matrices
M_k = [1; 1; 1];  % Uncertainty mapping
N_Fk = [0, 1, 1]; % Uncertainty distribution in F
N_Gk = [1, 1];    % Uncertainty distribution in G

%% Control Configuration

% Covariance matrices
[n_size, m_size] = size(G);
Q = NUSP * eye(n_size); % State weighting matrix
R = 0.01 * eye(m_size); % Control weighting matrix
P = eye(n_size);        % Ricatti initial matrix for LQR
P_r = eye(n_size);      % Ricatti initial matrix for RLQR

% Simulation steps
N = 50; 

% State vector
X_nom_lqr = zeros(n_size, N+1);
X_unc_lqr = zeros(n_size, N+1);
X_unc_rlqr = zeros(n_size, N+1);

% Control input
u_nom_lqr = zeros(m_size, N);
u_unc_lqr = zeros(m_size, N);
u_unc_rlqr = zeros(m_size, N);

% Initial conditions
Initial_State = [1.345; 1.75; 1.45];
X_nom_lqr(:,1) = Initial_State;
X_unc_lqr(:,1) = Initial_State;
X_unc_rlqr(:,1) = Initial_State;

%% Simulation - Item 1.1)

% Uncertainty balance and application
Delta_i = 0.75*(2 * rand - 1); % Uncertainty weight
DeltaF = M_k * Delta_i * N_Fk; % Uncertainty added to matrix F
DeltaG = M_k * Delta_i * N_Gk; % Uncertainty added to matrix G

% LQR Backward calculation
for k = 1:N+1
    [L, K, P] = standard_lqr(F, G, Q, R, P);
    [L_r, K_r, P_r] = robust_lqr(F, G, N_Fk, N_Gk, M_k, Q, R, P_r);
end

% Simulation loop
for k = 1:N
    u_nom_lqr(:,k) = K * X_nom_lqr(:,k);
    u_unc_lqr(:,k) = K * X_unc_lqr(:,k);
    u_unc_rlqr(:,k) = K_r * X_unc_rlqr(:,k);

    X_nom_lqr(:,k+1) = F * X_nom_lqr(:,k) + G * u_nom_lqr(:,k);
    X_unc_lqr(:,k+1) = (F + DeltaF) * X_unc_lqr(:,k) + (G + DeltaG) * u_unc_lqr(:,k);
    X_unc_rlqr(:,k+1) = (F + DeltaF) * X_unc_rlqr(:,k) + (G + DeltaG) * u_unc_rlqr(:,k);
end

%% Figure configuration

% Figure position and size
Figure_config = struct('PaperPositionMode', 'auto', ...
                      'Units', 'centimeters', ...
                      'Position', [9 9 25 15]);

% Font size and linewidth
LabelFontSizeValue = 12;
AxeslFontSizeValue = 8;
Object_config = struct('LineWidth', 2.50);

%% Output

% Ricatti P and Gain K Matrices Converged - Item 1.2)
disp("Matrix P Converged:");
disp(P_r)
disp("Matrix K Converged:");
disp(K_r)

% States comparison - Item 1.3)
% X_1
x_fig = figure();
set(x_fig, Figure_config);
subplot(3,1,1);
x1_obj = plot(1:N+1, X_nom_lqr(1,:), '--', ...
              1:N+1, X_unc_lqr(1,:), ...
              1:N+1, X_unc_rlqr(1,:));
set(x1_obj, Object_config);
legend('Nominal (LQR)', 'Uncertain (LQR)','Uncertain (RLQR)')
ylabel('$x_1$', 'Interpreter', 'Latex', 'fontsize', LabelFontSizeValue)
axis tight
ax = gca;
ax.YGrid = 'on';
ax.FontSize = LabelFontSizeValue;
title('States comparison')

% X_2
subplot(3,1,2);
x2_obj = plot(1:N+1, X_nom_lqr(2,:), '--', ...
              1:N+1, X_unc_lqr(2,:), ...
              1:N+1, X_unc_rlqr(2,:));
set(x2_obj, Object_config);
legend('Nominal (LQR)', 'Uncertain (LQR)','Uncertain (RLQR)')
ylabel('$x_2$', 'Interpreter', 'Latex', 'fontsize', LabelFontSizeValue)
axis tight
ax = gca;
ax.YGrid = 'on';
ax.FontSize = LabelFontSizeValue;

% X_3
subplot(3,1,3);
x3_obj = plot(1:N+1, X_nom_lqr(3,:), '--', ...
              1:N+1, X_unc_lqr(3,:), ...
              1:N+1, X_unc_rlqr(3,:));
set(x3_obj, Object_config);
legend('Nominal (LQR)', 'Uncertain (LQR)','Uncertain (RLQR)')
ylabel('$x_3$', 'Interpreter', 'Latex', 'fontsize', LabelFontSizeValue)
xlabel('k', 'Interpreter', 'Latex', 'fontsize', LabelFontSizeValue)
axis tight
ax = gca;
ax.YGrid = 'on';
ax.FontSize = LabelFontSizeValue;

LabelFontSizeValue = LabelFontSizeValue * 2;
AxeslFontSizeValue = AxeslFontSizeValue * 2;

% Control inputs - Item 1.4)
% u_1
Control_fig = figure();
set(Control_fig, Figure_config);
subplot(2,1,1);
u1_plot_obj = plot(1:N, u_nom_lqr(1,:), '--', ...
                   1:N, u_unc_lqr(1,:), ...
                   1:N, u_unc_rlqr(1,:));
set(u1_plot_obj, Object_config);
title('Control input')
ylabel('$u_1$', 'Interpreter', 'Latex', 'fontsize', LabelFontSizeValue)
axis tight
ax = gca;
ax.YGrid = 'on';
ax.FontSize = AxeslFontSizeValue;
legend('Nominal (LQR)', 'Uncertain (LQR)', 'Uncertain (RLQR)', 'Location', 'Southeast')

% u_2
subplot(2,1,2);
u2_plot_ob = plot(1:N, u_nom_lqr(2,:), '--', ...
                  1:N, u_unc_lqr(2,:), ...
                  1:N, u_unc_rlqr(2,:));
set(u2_plot_ob, Object_config);
ylabel('$u_2$', 'Interpreter', 'Latex', 'fontsize', LabelFontSizeValue)
xlabel('k', 'Interpreter', 'Latex', 'fontsize', LabelFontSizeValue)
axis tight
ax = gca;
ax.YGrid = 'on';
ax.FontSize = AxeslFontSizeValue;
legend('Nominal (LQR)', 'Uncertain (LQR)', 'Uncertain (RLQR)', 'Location', 'Southeast')

% Eigenvalues Plot - Item 1.5)
Openloop_values = eig((F+DeltaF));
Closedloop_values = eig((F+DeltaF)+(G+DeltaG)*K_r);

Eigenvalues_fig = figure;
set(Eigenvalues_fig, Figure_config);
th = 0:pi/100:2*pi;
plot(cos(th), sin(th), 'LineWidth', 2);
hold on;
plot(real(Openloop_values), imag(Openloop_values), '.', ...
    'MarkerSize', 25, 'LineWidth', 2);
plot(real(Closedloop_values), imag(Closedloop_values), '.', ...
    'MarkerSize', 25, 'LineWidth', 2);
grid on;
axis equal;
legend('', 'Open-loop', 'Closed-loop (RLQR)','Location', 'Southeast')
title('Eigenvalues')
