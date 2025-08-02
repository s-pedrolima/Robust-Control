%% Initialization
clear vars; clc; close all;

% Truck parameters (Based on Table 4 and system_data from PDF)
a1 = 1.734;     b1 = 2.415;
a2 = 4.8;       b2 = 3.2;
l1 = a1 + b1;   l2 = a2 + b2; % l1 and l2 from parameters
h1 = b1 - 0.29; % Based on h1 = b1 + et1 where et1 = -0.29
v = 16.667;     m1 = 8909;
payload = 24000; m2 = 9370 + payload;
J1 = 41566;     J2 = 404360;
c1 = 345155;    c2 = 927126;    c3 = 1158008;
g = 9.8;

% Disturbance (from your original code)
a = 10; b = 0.1;

% System matrices (Eq. 38 from your problem description)
M = [m1+m2,      -m2*(h1+a2),    -m2*a2,         0, 0, 0;
     -m2*h1,      J1+m2*h1*(h1+a2), m2*h1*a2,     0, 0, 0;
     -m2*a2,      m2*a2*(h1+a2),   J2+m2*a2^2,    0, 0, 0;
     0, 0, 0, 1, 0, 0;
     0, 0, 0, 0, 1, 0;
     0, 0, 0, 0, 0, 1];

A1 = [-(c1+c2+c3)/v, (c3*(h1+l2)-a1*c1+b1*c2-(m1+m2)*v^2)/v, c3*l2/v, c3, 0, 0;
      (c3*h1-a1*c1+b1*c2)/v, (m2*h1*v^2-a1^2*c1-b1^2*c2-c3*h1*(h1+l2))/v, -c3*h1*l2/v, -c3*h1, 0, 0;
      c3*l2/v, (m2*a2*v^2-c3*l2*(h1+l2))/v, -c3*l2^2/v, -c3*l2, 0, 0;
      0, 0, 1, 0, 0, 0;
      1, 0, 0, 0, 0, v;
      0, 1, 0, 0, 0, 0];
      
B1 = [c1; a1*c1; 0; 0; 0; 0];
Bw = [1; 1; 1; 0; 0; 0];  % Disturbance matrix

% Discretization (Tustin)
Ts = 0.01;
sysc = ss(M\A1, M\[B1, Bw], eye(6), 0);
sysd = c2d(sysc, Ts, 'tustin');
F = sysd.A;
G = sysd.B(:,1);   % Control matrix
D = sysd.B(:,2);   % Disturbance matrix

%% Reference Trajectory Generation
Tf = 30; t = 0:Ts:Tf; N = length(t);
u_ref = zeros(1, N);
t_lc = 0:Ts:5; % 5-second maneuver time

% Define reference input for a double lane change (similar to PDF)
start_1 = 10 / Ts; % Start at 10s
end_1 = start_1 + length(t_lc) -1;
u_ref(start_1:end_1) = 0.01*sin(0.4*pi*t_lc);

start_2 = 20 / Ts; % Start at 20s
end_2 = start_2 + length(t_lc) -1;
u_ref(start_2:end_2) = -0.01*sin(0.4*pi*t_lc);

% Generate reference states
x_ref = zeros(size(F,1), N);
for k = 1:N-1
    x_ref(:,k+1) = F*x_ref(:,k) + G*u_ref(k);
end

% Generate global reference path
[posex_ref, posey_ref] = get_global_path(x_ref, v, Ts);


%% Controller Design
% Common settings
n = size(F,1); m = size(G,2);
P0 = eye(n);
max_iter = 1000;
tol = 1e-6;

% Replace S with the sum of the digits of your USP
S = 1+0+7+4+8+1+9+8;  

% --- Início da Modificação ---

fprintf('Iniciando o projeto dos controladores...\n');
fprintf('-----------------------------------------\n');

% Case 1: Q1 = 0.1I
[K_hinf1, gamma1] = design_hinf(F, G, D, 0.1*eye(n), 0.1, 10);
fprintf('Gamma mínimo para Q = 0.1*I: %.4f\n', gamma1);

% Case 2: Q2 = SI
[K_hinf2, gamma2] = design_hinf(F, G, D, S*eye(n), 0.1, 10);
fprintf('Gamma mínimo para Q = S*I (S=%d): %.4f\n', S, gamma2);

% Case 3: Q3 = 10000I
[K_hinf3, gamma3] = design_hinf(F, G, D, 10000*eye(n), 1, 100);
fprintf('Gamma mínimo para Q = 10000*I: %.4f\n', gamma3);

% Case 4: Fixed gamma = 1000
[~, K_hinf4] = solve_hinf_riccati(F, G, D, S*eye(n), 1000, max_iter, tol);
fprintf('Gamma fixo para Q = S*I: 1000 (por definição)\n');

% LQR (case 2 as reference)
[K_lqr, ~] = dlqr(F, G, S*eye(n), 1);
fprintf('Controlador LQR projetado.\n');
fprintf('-----------------------------------------\n');

%% Simulation
x0 = [0; 0; 0; 0; 0.3; -0.1];  % Initial condition

% Disturbance generation
w = a * exp(-b*t) .* cos(2*t);

% Simulation for each controller
controllers = {K_hinf1, K_hinf2, K_hinf3, K_hinf4, K_lqr};
labels = {'H-inf Q=0.1I', 'H-inf Q=SI', 'H-inf Q=10000I', 'H-inf gamma=1000', 'LQR'};
results = cell(1,5);

for ctrl = 1:length(controllers)
    K = controllers{ctrl};
    x = zeros(n,N); u = zeros(m,N-1);
    e = zeros(n,N);
    x(:,1) = x0;
    
    for k = 1:N-1
        % Path-following control law
        e(:,k) = x(:,k) - x_ref(:,k);
        u(:,k) = -K * e(:,k) + u_ref(k);
        
        % Update state with disturbance
        x(:,k+1) = F*x(:,k) + G*u(:,k) + D*w(k);
    end
    
    % Calculate global path for the simulated vehicle
    [posex, posey] = get_global_path(x, v, Ts);
    results{ctrl} = struct('x',x, 'u',u, 't',t, 'posex', posex, 'posey', posey, 'e', e);
end

%% Plotting Results
plot_styles = {'-r', '-g', '-b', '-c', '--m'};
state_labels = {'Lateral Velocity y_1_dot (m/s)', 'Yaw Rate psi_dot (rad/s)', ...
                'Articulation Rate phi_dot (rad/s)', 'Articulation Angle phi (rad)', ...
                'Displacement Error rho (m)', 'Orientation Error theta (rad)'};

% --- Plot 1: Global Position (Já estava correto)
figure('Name', 'Global Position', 'Position', [100 100 900 600]);
hold on;
plot(posex_ref, posey_ref, 'k:', 'LineWidth', 2.5, 'DisplayName', 'Reference');
for ctrl = 1:length(controllers)
    plot(results{ctrl}.posex, results{ctrl}.posey, plot_styles{ctrl}, 'LineWidth', 1.5, 'DisplayName', labels{ctrl});
end
hold off;
title('Global Position');
xlabel('X (m)');
ylabel('Y (m)');
legend('show', 'Location', 'northwest');
grid on;
xlim([0 500]);

% --- Plot 2: Control Input (Já estava correto)
figure('Name', 'Control Input', 'Position', [100 100 900 400]);
hold on;
plot(t, u_ref, 'k:', 'LineWidth', 2.5, 'DisplayName', 'u_{ref}');
for ctrl = 1:length(controllers)
    plot(results{ctrl}.t(1:end-1), results{ctrl}.u, plot_styles{ctrl}, 'LineWidth', 1.5, 'DisplayName', labels{ctrl});
end
hold off;
title('Control Input \alpha (rad)');
xlabel('Time (s)');
ylabel('\alpha (rad)');
legend('show');
grid on;

% --- Plots 3, 4, 5, 6: Estados 1 a 4
for i = [1, 2, 3, 4] % Plota apenas os 4 primeiros estados
    figure('Name', state_labels{i}, 'Position', [100 100 900 400]);
    hold on;
    plot(t, x_ref(i,:), 'k:', 'LineWidth', 2.5, 'DisplayName', 'Reference');
    for ctrl = 1:length(controllers)
        plot(results{ctrl}.t, results{ctrl}.x(i,:), plot_styles{ctrl}, 'LineWidth', 1.5, 'DisplayName', labels{ctrl});
    end
    hold off;
    title(state_labels{i});
    xlabel('Time (s)');
    ylabel(state_labels{i});
    legend('show');
    grid on;
end

% --- Plot 7: Gráfico de ERRO do estado 5 (rho)
figure('Name', state_labels{5}, 'Position', [100 100 900 400]);
hold on;
% A referência para o ERRO é zero
plot(t, zeros(size(t)), 'k:', 'LineWidth', 2.5, 'DisplayName', 'Reference (Error=0)');
for ctrl = 1:length(controllers)
    % Plota o 5º elemento do vetor de erro 'e'
    plot(results{ctrl}.t, results{ctrl}.e(5,:), plot_styles{ctrl}, 'LineWidth', 1.5, 'DisplayName', labels{ctrl});
end
hold off;
title(state_labels{5});
xlabel('Time (s)');
ylabel(state_labels{5});
legend('show');
grid on;

% --- Plot 8: Estado 6 (theta)
i = 6;
figure('Name', state_labels{i}, 'Position', [100 100 900 400]);
hold on;
plot(t, x_ref(i,:), 'k:', 'LineWidth', 2.5, 'DisplayName', 'Reference');
for ctrl = 1:length(controllers)
    plot(results{ctrl}.t, results{ctrl}.x(i,:), plot_styles{ctrl}, 'LineWidth', 1.5, 'DisplayName', labels{ctrl});
end
hold off;
title(state_labels{i});
xlabel('Time (s)');
ylabel(state_labels{i});
legend('show');
grid on;


%% Helper Functions

function [posex, posey] = get_global_path(x_data, v, Ts)
    % Para manobras de desvio de pista com velocidade constante 'v',
    % podemos aproximar a posição X pela integração da velocidade.
    % O estado 'rho' (5º estado) representa a posição lateral Y.
    N = size(x_data, 2);
    t = (0:N-1) * Ts;

    posex = v * t;       % Posição X é a velocidade constante * tempo
    posey = x_data(5,:); % Posição Y é o próprio estado rho
end

function [K_inf, gamma_opt] = design_hinf(F, G, D, Q, gamma_min, gamma_max)
    % Finds smallest viable gamma
    gamma_test = linspace(gamma_min, gamma_max, 50);
    viable = false;
    
    for g = gamma_test
        try
            [~, K] = solve_hinf_riccati(F, G, D, Q, g, 1000, 1e-6);
            gamma_opt = g;
            K_inf = K;
            viable = true;
            break;
        catch
            continue;
        end
    end
    
    if ~viable
        error('No viable γ found in the range');
    end
end

function [P_inf, K_inf] = solve_hinf_riccati(F, G, D, Q, gamma, max_iter, tol)
    n = size(F,1);
    P_prev = eye(n);
    I = eye(n);
    
    for i = 1:max_iter
        [Lambda, P_next] = h_inf_ga(F, G, D, Q, P_prev, gamma);
        if norm(P_next - P_prev, 'fro') < tol
            break;
        end
        P_prev = P_next;
    end
    
    P_inf = P_next;
    % K_inf = -(G' * P_inf / Lambda * F);  % Equivalent to inv(Lambda)
    K_inf = inv(G'*P_inf*G + 1) * G' * P_inf * F; % More stable formulation
end

function [Lambda, P] = h_inf_ga(F, G, D, Q, Pp, gamma)
    n = size(Pp,1);
    I = eye(n);
    
    disturbance_term = (1/gamma^2) * (D * D');
    control_term = G * G';
    
    % Calculate P_k+1 using the Riccati equation
    Term_inv = inv(I + (control_term - disturbance_term)*Pp);
    P = Q + F' * Pp * Term_inv * F;
    
    % Stability check
    Xi = gamma^2 * eye(size(D,2)) - D' * P * D;
    if any(eig(Xi) <= 0)
        error('Gamma too small: Riccati condition violated');
    end
    
    % Lambda is not explicitly needed for K calculation in this form
    Lambda = []; 
end