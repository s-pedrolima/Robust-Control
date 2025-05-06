% Controle de posição de um pêndulo invertido sobre um carrinho
% Seguindo uma trajetória pré-determinada
% Considerando a variação contínua da massa ao longo do tempo

close all;
clearvars;
clc;
import Controllers.*
import Aux.*

%% Parâmetros e Configurações do Sistema

g = 9.81; % Aceleração gravitacional
I_n = 0.099; % Momento de inércia
M0 = 2.4;   % Massa inicial do carrinho
m0 = 0.23;  % Massa inicial do pêndulo
l = 0.4; % Comprimento da haste
b = 0.05; % Coeficiente de atrito
Xi = 0.005; % Coeficiente de amortecimento
k_m = 8; % Constante do motor
Ts = 10^-2; % Tempo de amostragem: 10 ms

% Número de amostras
N = 50/Ts;

% Matrizes de covariância
n_size = 4; % Dimensão do estado
m_size = 1; % Dimensão do controle
P = eye(n_size); % Matriz de Ricatti inicial
Q = eye(n_size); % Penalização dos estados
R = eye(m_size); % Penalização de controle

% Inicialização das variáveis de estado
X = zeros(n_size, N+1);
u = zeros(m_size, N);

% Posição inicial da haste
xa = 0; % Posição linear do carrinho
theta = 0; % Posição angular da haste
X(:,1) = [theta; xa; 0; 0]; 

% Referência de posição - Caminho a ser seguido
Amplitude = 2 * (0.35);
Signal_period = 30;
ref_signal = Amplitude * (-gensig('square', Signal_period, N*Ts, Ts) + 0.5);

% Vetor de tempo
time = (0:N) * Ts;

%% Loop Principal - Aplicação de Controle
for k = 1:N
    t = (k-1) * Ts;  % Tempo atual para cálculo das massas
    
    % Atualização da massa do carrinho e do pêndulo
    M_t = m0 * (1 + 0.1*sin(0.2*t));
    m_t = M0 * (1 + 0.05*cos(0.1*t));
    
    % Atualização das matrizes do sistema para o tempo atual
    [F, G, H, D] = generate_pendulum_model(g, I_n, M_t, m_t, l, b, Xi, k_m, Ts);
    
    % Cálculo do LQR
    [L, K, P] = standard_lqr(F, G, Q, R, P);
   
    % Entrada de controle ótimo
    u(:,k) = K * X(:,k);

    % Saturação do controle (máx. ±50 N)
    u(:,k) = max(min(u(:,k), 50), -50);

    % Caminho a ser seguido
    Position_ref = ref_signal(k) * [0; 0; 0; Ts];

    % Evolução do sistema
    X(:,k+1) = F * X(:,k) + G * u(:,k) + Position_ref;

    % Limitar o carrinho ao trilho de 2 metros (±1 m)
    if abs(X(2,k+1)) > 1
        fprintf('🚨 Carrinho saiu dos trilhos no tempo %.2f s\n', t+Ts);
        break;
    end
    
    % Se a haste caiu, ou seja, se theta (assumindo theta = X(1,k+1)) ultrapassar um limiar
    if abs(X(1,k+1)) > pi
        fprintf('🚨 Haste caiu no tempo %.2f s\n', t+Ts);
        break;
    end
end

%% Gráficos - Ajustando tamanho de fontes e espessura das curvas
figure;

fontSize = 14; % Tamanho da fonte
lineWidth = 2; % Espessura da linha

% Posição Linear
subplot(3,1,1);
plot(time, ref_signal, 'k--', ...
     time, X(2,:), 'r-', 'LineWidth', lineWidth);
title('Linear Position', 'FontSize', fontSize+2, 'FontWeight', 'bold');
ylabel('Position (m)', 'FontSize', fontSize);
legend('Reference', 'Car Path', 'Location', 'best');
grid on;
xlim([0 max(time)]);
set(gca, 'FontSize', fontSize);

% Posição Angular
subplot(3,1,2);
plot(time, X(1,:), 'b', 'LineWidth', lineWidth);
title('Angular Position', 'FontSize', fontSize+2, 'FontWeight', 'bold');
ylabel('Angle (rad)', 'FontSize', fontSize);
grid on;
xlim([0 max(time)]);
set(gca, 'FontSize', fontSize);

% Sinal de Controle
subplot(3,1,3);
plot(time(1:end-1), u(1,:), 'g', 'LineWidth', lineWidth);
title('Control Input', 'FontSize', fontSize+2, 'FontWeight', 'bold');
ylabel('Input (N)', 'FontSize', fontSize);
xlabel('Time (s)', 'FontSize', fontSize);
grid on;
xlim([0 max(time)]);
set(gca, 'FontSize', fontSize);