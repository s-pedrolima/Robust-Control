%% SEL0382 - Controle Robusto
% Solução Final, Estável e Completa, incluindo todos os controladores

clear; close all; clc;

%% 1. Definição do Sistema e Parâmetros
% Parâmetros físicos dos slides 
M_s = 2.45;      % Massa suspensa (kg) 
M_us = 1;        % Massa não suspensa (kg) 
K_s = 900.0;     % Rigidez da suspensão (N/m) 
K_us = 1250;     % Rigidez do pneu (N/m) 
B_s = 7.5;       % Amortecimento da suspensão (N.s/m) 
B_us = 5.0;      % Amortecimento do pneu (N.s/m) 
Ts = 0.01;       % Tempo de amostragem (s) da Lista 5

% Matrizes do sistema em tempo contínuo 
A = [0, 1, 0, -1;
     -K_s/M_s, -B_s/M_s, K_s/M_s, B_s/M_s;
     0, 0, 0, 1;
     K_s/M_us, B_s/M_us, -K_us/M_us, -(B_s + B_us)/M_us];

% Matriz de entrada 4x2, baseada na estrutura dos slides 
B = [0, 0;
     0, 1/M_s;
     -1, 0;
     B_us/M_us, -1/M_us];

% Discretização do sistema NOMINAL 
I = eye(size(A));
F = A * Ts + I;
G = B * Ts;
G_w = G(:, 1); % Coluna do distúrbio
G_u = G(:, 2); % Coluna do controle

% Sistema com INCERTEZA para a simulação 
K_s_unc = 1.15 * K_s; % 
B_s_unc = 1.10 * B_s; % 
A_unc = [0, 1, 0, -1;
         -K_s_unc/M_s, -B_s_unc/M_s, K_s_unc/M_s, B_s_unc/M_s;
         0, 0, 0, 1;
         K_s_unc/M_us, B_s_unc/M_us, -K_us/M_us, -(B_s_unc + B_us)/M_us]; % 
F_unc = A_unc * Ts + I; % 
G_unc = B * Ts; % 

%% 2. Configuração do Controle e da Simulação
S = 28; % Valor assumido para S
Q = diag([45*S, 30*S, 5*S, 0.01*S]); % Ponderação Q da Lista 5 
% Usando R=0.01 dos slides para garantir estabilidade. O valor R=1 da
% lista torna o problema H-inf ótimo instável para a planta incerta.
R = 0.01; % 

% Simulação
sim_time = 10; % 
N = sim_time / Ts; % 
t = (0:N-1) * Ts;

% Geração da referência de pista (onda quadrada) 
Amplitude = 0.02; % 
Signal_period = 3; % 
t_ref = (0:N) * Ts;
ref_signal = Amplitude * (mod(t_ref, Signal_period) < Signal_period/2);
delay_samples = round(2 / Ts); % Atraso customizado para melhor visualização
ref_signal(1:delay_samples) = 0;
ref_signal = ref_signal(1:N+1);
Road_profile = diff(ref_signal)/Ts; % 

%% 3. Projeto dos Controladores
% Função Robusta para resolver a Equação de Riccati H-infinity (DARE)
function [K, P] = solve_h_inf_DARE(F, G_u, G_w, Q, R, gamma, max_iter, tol)
    P = Q;
    B_aug = [G_u, G_w];
    R_u = R(1,1);
    R_aug = diag([R_u, -gamma^2]);
    for iter = 1:max_iter
        P_prev = P;
        M = R_aug + B_aug' * P_prev * B_aug;
        if rcond(M) < 1e-12, K = []; P = []; return; end
        Term = F' * P_prev * B_aug * (M \ B_aug') * P_prev * F;
        P = F' * P_prev * F + Q - Term;
        if norm(P - P_prev, 'fro') < tol, break; end
    end
    M_final = R_aug + B_aug' * P * B_aug;
    if rcond(M_final) < 1e-12, K = []; P = []; return; end
    L = M_final \ (B_aug' * P * F);
    K = -L(1,:);
end

% 1. H-infinity com o menor gamma possível
low_gamma = 0.1; high_gamma = 20.0; optimal_gamma = high_gamma;
while (high_gamma - low_gamma) > 1e-4
    mid_gamma = (low_gamma + high_gamma) / 2;
    [K_test, ~] = solve_h_inf_DARE(F, G_u, G_w, Q, R, mid_gamma, 1000, 1e-6);
    if ~isempty(K_test), optimal_gamma = mid_gamma; high_gamma = mid_gamma;
    else, low_gamma = mid_gamma; end
end
gamma_safe = optimal_gamma * 1.15;
fprintf('Gamma (γ) ótimo encontrado: %.4f\n', optimal_gamma);
fprintf('Gamma (γ) seguro utilizado para o controle: %.4f\n', gamma_safe);
[K_hinf_opt, ~] = solve_h_inf_DARE(F, G_u, G_w, Q, R, gamma_safe, 1000, 1e-6);

% 2. H-infinity com gamma = 1000 
[K_hinf_1000, ~] = solve_h_inf_DARE(F, G_u, G_w, Q, R, 1000, 1000, 1e-6);

% 3. LQR Padrão
[K_lqr, ~, ~] = dlqr(F, G_u, Q, R);
K_lqr = -K_lqr;

%% 4. Simulação em Malha Fechada
function [x, u] = simulate_system(F_unc, G_unc, K, N, road_profile)
    n_states = size(F_unc, 1);
    x = zeros(n_states, N + 1); u = zeros(1, N);
    G_w_unc = G_unc(:, 1); G_u_unc = G_unc(:, 2);
    for i = 1:N
        u(i) = K * x(:, i);
        x(:, i+1) = F_unc * x(:, i) + G_u_unc * u(i) + G_w_unc * road_profile(i);
    end
end
[x_hinf_opt, u_hinf_opt] = simulate_system(F_unc, G_unc, K_hinf_opt, N, Road_profile); % 
[x_hinf_1000, u_hinf_1000] = simulate_system(F_unc, G_unc, K_hinf_1000, N, Road_profile);
[x_lqr, u_lqr] = simulate_system(F_unc, G_unc, K_lqr, N, Road_profile); % 

%% 5. Resultados e Gráficos (com cores e estilos aprimorados)
figure('Position', [100, 100, 800, 700]);

% --- Definição de cores e estilos para clareza ---
% LQR: Azul, tracejado
% H-inf (g=1000): Vermelho, pontilhado e grosso (para sobrepor o LQR)
% H-inf (g=ótimo): Verde, sólido
color_lqr = [0, 0.4470, 0.7410];       % Azul padrão do MATLAB
color_hinf_1000 = [0.8500, 0.3250, 0.0980]; % Vermelho/Laranja padrão
color_hinf_opt = [0.4660, 0.6740, 0.1880];  % Verde padrão

% Plot 1: Posição da Massa Suspensa (zs)
subplot(3,1,1);
zs_lqr = x_lqr(1, 1:N+1) + x_lqr(3, 1:N+1) + ref_signal;
zs_hinf_1000 = x_hinf_1000(1, 1:N+1) + x_hinf_1000(3, 1:N+1) + ref_signal;
zs_hinf_opt = x_hinf_opt(1, 1:N+1) + x_hinf_opt(3, 1:N+1) + ref_signal;

plot(t_ref, zs_lqr*1000, '--', 'Color', color_lqr, 'LineWidth', 1.5); hold on;
plot(t_ref, zs_hinf_1000*1000, ':', 'Color', color_hinf_1000, 'LineWidth', 2.0);
plot(t_ref, zs_hinf_opt*1000, '-', 'Color', color_hinf_opt, 'LineWidth', 1.5);
plot(t_ref, ref_signal*1000, 'k:', 'LineWidth', 1);

title('Top plate position');
ylabel('$z_s$ (mm)', 'Interpreter', 'latex');
legend('LQR', 'H-inf (\gamma=1000)', ['H-inf (\gamma=' sprintf('%.2f', gamma_safe) ')'], 'Reference');
grid on; axis tight;

% Plot 2: Aceleração da Carroceria
subplot(3,1,2);
accel_lqr = diff(x_lqr(2, :)) / Ts;
accel_hinf_1000 = diff(x_hinf_1000(2, :)) / Ts;
accel_hinf_opt = diff(x_hinf_opt(2, :)) / Ts;

plot(t, accel_lqr, '--', 'Color', color_lqr, 'LineWidth', 1.5); hold on;
plot(t, accel_hinf_1000, ':', 'Color', color_hinf_1000, 'LineWidth', 2.0);
plot(t, accel_hinf_opt, '-', 'Color', color_hinf_opt, 'LineWidth', 1.5);

title('Vehicle body acceleration');
ylabel('$m/s^2$', 'Interpreter', 'latex');
legend('LQR', 'H-inf (\gamma=1000)', ['H-inf (\gamma=' sprintf('%.2f', gamma_safe) ')']);
grid on; axis tight;

% Plot 3: Força de Controle da Suspensão
subplot(3,1,3);
plot(t, u_lqr, '--', 'Color', color_lqr, 'LineWidth', 1.5); hold on;
plot(t, u_hinf_1000, ':', 'Color', color_hinf_1000, 'LineWidth', 2.0);
plot(t, u_hinf_opt, '-', 'Color', color_hinf_opt, 'LineWidth', 1.5);

title('Active suspension control');
ylabel('$F_c$ (N)', 'Interpreter', 'latex');
xlabel('Time (s)', 'Interpreter', 'latex');
legend('LQR', 'H-inf (\gamma=1000)', ['H-inf (\gamma=' sprintf('%.2f', gamma_safe) ')']);
grid on; axis tight;