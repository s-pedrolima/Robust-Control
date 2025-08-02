%% H-infinity control of a Cart-pendulum Example

% Initialization
close all;
clearvars;
clc;
import Controllers.*

%% System Parameters
Ts = 0.0001;

% System matrices
F = [1 -4.9*10^-6 2.5*10^-3 0;
     -1.2*10^-6 1 1.6*10^-4 0;
     9.9*10^-5 0 1 0;
     0 -9.9*10^-5 0 1];
     
G = [-7.82*10^-4;
      3.53*10^-4;
      0;
      0];
  
D = Ts*[1;
        1;
        0;
        0];

%% Weighting Matrices
[n_size, m_size] = size(G);
P_hinf = 1e3 * eye(n_size);
P_lqr = 1e3 * eye(n_size);
Q = 1e2 * eye(n_size);
R = 1e-2 * eye(m_size);

%% Simulation Setup
N = 60/Ts; % Number of steps

% State vectors
X_lqr = zeros(n_size, N+1);
X_hinf = zeros(length([6; 1e5]), n_size, N+1);

% Control inputs
u_lqr = zeros(m_size, N);
u_hinf = zeros(length([6; 1e5]), m_size, N);

% Cost vectors
Cost_LQR = zeros(1, N);
Cost_Hinf = zeros(length([6; 1e5]), N);

% Reference signal
Amplitude = 2 * (0.35);
Signal_period = 30;
ref_signal = Amplitude * (-gensig('square', Signal_period, N*Ts, Ts) + 0.5);

%% Controller Design
% Disturbance
w = D*randn(1, N);

% Gamma values for H_inf
gamma_values = [6; 1e5];
gamma_length = length(gamma_values);
P_hinf_hist = zeros(gamma_length, n_size, n_size);
legend_vector = cell(gamma_length + 1, 1); % +1 for LQR

% LQR calculation
for i = 1:N
    [~, K, P_lqr] = standard_lqr(F, G, Q, R, P_lqr);
end

% H-infinity controller calculation
for k = 1:gamma_length
    gamma = gamma_values(k);
    
    for i = 1:N
        [Lambda, P_hinf] = h_inf_ga(F, G, D, Q, P_hinf, gamma);
    end
    
    legend_vector{k} = sprintf('$\\gamma = %s$', num2str(gamma));
    P_hinf_hist(k, :, :) = P_hinf;
    P_hinf = 1e3 * eye(n_size); % Reset P for next iteration
end

%% Simulation
for i = 1:N
    Position_ref = ref_signal(i) * [0; 0; 0; Ts];
    
    % H-infinity control
    for k = 1:gamma_length
        u_hinf(k, :, i) = -G' * squeeze(P_hinf_hist(k, :, :)) * inv(Lambda) * F * squeeze(X_hinf(k, :, i))';
        X_hinf(k, :, i+1) = F * squeeze(X_hinf(k, :, i))' + G * squeeze(u_hinf(k, :, i)) + w(:, i) + Position_ref;
        Cost_Hinf(k, i) = squeeze(X_hinf(k, :, i)) * squeeze(P_hinf_hist(k, :, :)) * squeeze(X_hinf(k, :, i))';
    end
    
    % LQR control
    u_lqr(:, i) = K * X_lqr(:, i);
    X_lqr(:, i+1) = F * X_lqr(:, i) + G * u_lqr(:, i) + w(:, i) + Position_ref;
    Cost_LQR(i) = X_lqr(:, i)' * P_lqr * X_lqr(:, i);
end

legend_vector{end} = sprintf('$LQR$');

%% Figure Configuration
Figure_config = struct('PaperPositionMode', 'auto', ...
                      'Units', 'centimeters', ...
                      'Position', [8 4 40 22]);
                  
LabelFontSizeValue = 20;
AxeslFontSizeValue = 16;
Object_config = struct('LineWidth', 2.50);

%% Plotting Results
% State Plots
states_fig = figure(1);
set(states_fig, Figure_config);

% Angular position
subplot(3, 1, 1)
x1_obj = plot((1:N+1)*Ts, squeeze(X_hinf(:, 1, :)), ...
             (1:N+1)*Ts, X_lqr(1, :), '-.');
set(x1_obj, Object_config);
ylabel('$rad$', 'Interpreter', 'Latex', 'fontsize', LabelFontSizeValue)
xlabel('t (s)', 'Interpreter', 'Latex', 'fontsize', LabelFontSizeValue)
title('Angular position')
grid on
ax = gca;
ax.FontSize = AxeslFontSizeValue;
axis tight

% Linear position
subplot(3, 1, 2)
x2_obj = plot((1:N+1)*Ts, squeeze(X_hinf(:, 2, :)), ...
             (1:N+1)*Ts, X_lqr(2, :), '-.', ...
             (1:N+1)*Ts, ref_signal, 'k:');
set(x2_obj, Object_config);
ylabel('$m$', 'Interpreter', 'Latex', 'fontsize', LabelFontSizeValue)
xlabel('t (s)', 'Interpreter', 'Latex', 'fontsize', LabelFontSizeValue)
legend(legend_vector, 'Interpreter', 'Latex', 'fontsize', LabelFontSizeValue)
title('Linear position')
grid on
ax = gca;
ax.FontSize = AxeslFontSizeValue;
axis tight

% Control input
subplot(3, 1, 3)
u1_plot_obj = plot((1:N)*Ts, squeeze(u_hinf(:, 1, :)), ...
                  (1:N)*Ts, u_lqr(1, :), '-.');
set(u1_plot_obj, Object_config);
ylabel('$u$', 'Interpreter', 'Latex', 'fontsize', LabelFontSizeValue)
xlabel('t (s)', 'Interpreter', 'Latex', 'fontsize', LabelFontSizeValue)
title('Control input')
grid on
ax = gca;
ax.FontSize = AxeslFontSizeValue;
axis tight

% Cost Plot
Costs_fig = figure(2);
set(Costs_fig, Figure_config);
c_plot_obj = plot((1:N)*Ts, Cost_Hinf, ...
                 (1:N)*Ts, Cost_LQR, '-.');
set(c_plot_obj, Object_config);
title('Costs $H \infty$ and LQR', 'Interpreter', 'Latex', 'fontsize', LabelFontSizeValue)
xlabel('t (s)', 'Interpreter', 'Latex', 'fontsize', LabelFontSizeValue)
legend(legend_vector, 'Interpreter', 'Latex', 'fontsize', LabelFontSizeValue, 'Location', 'Southeast')
grid on
ax = gca;
ax.FontSize = AxeslFontSizeValue;
axis tight

% Individual Control Plots
for i = 1:gamma_length
    Control_fig = figure(2 + i);
    set(Control_fig, Figure_config);
    control_plot_obj = plot((1:N)*Ts, squeeze(u_hinf(i, 1, :)));
    set(control_plot_obj, Object_config);
    ylabel('$u$', 'Interpreter', 'Latex', 'fontsize', LabelFontSizeValue)
    xlabel('t (s)', 'Interpreter', 'Latex', 'fontsize', LabelFontSizeValue)
    legend(strcat('$\gamma = $', num2str(gamma_values(i))), 'Interpreter', 'Latex', 'fontsize', LabelFontSizeValue)
    xlim([0 60])
    ylim([-1 1])
    grid on
    ax = gca;
    ax.FontSize = AxeslFontSizeValue;
end

% LQR Control Plot
Control_fig = figure(2 + gamma_length + 1);
set(Control_fig, Figure_config);
control_plot_obj = plot((1:N)*Ts, u_lqr(1, :));
set(control_plot_obj, Object_config);
ylabel('$u$', 'Interpreter', 'Latex', 'fontsize', LabelFontSizeValue)
xlabel('t (s)', 'Interpreter', 'Latex', 'fontsize', LabelFontSizeValue)
legend('LQR', 'fontsize', LabelFontSizeValue)
xlim([0 60])
ylim([-1 1])
grid on
ax = gca;
ax.FontSize = AxeslFontSizeValue;