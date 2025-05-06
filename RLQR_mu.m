%% RLQR Numerical Example

clc;
clearvars;
close all; 
import Controllers.*
import Aux.*

%% System Configuration

% System matrices
F = [1.1, 0, 0;
     0, 0, 1.2;
    -1, 1, 0];

G = [0, 1;
     1, 1;
    -1, 0];

% Include unceartainties? Y(1)/N(0)
a = 1;

% Unceartainties matrices
H = [.7;
     .5;
    -.7 ];

E_F = a*[.4, .5, -.6];
E_G = a*[.4, -.4];

% Initial states
xo = [.1; -.1; .5];

% Penalty parameters: mu>0
mu = [1e2, 1e5, 1e10];

% Loop length
Horizonte = 70;

% Covariance matrices
[n_size, m_size] = size(G);
Q = eye(n_size, n_size);                % State weighting matrix
R = eye(m_size, m_size);                % Control weighting matrix
K = zeros(m_size, n_size, Horizonte);   % Gain matrix
L = zeros(n_size, n_size, Horizonte);   % L matrix
P = zeros(n_size, n_size, Horizonte+1); % Ricatti matrix
P(:, :, 1) = eye(n_size);               % Initial P value

% Figures settings
Titles = ["Poles: Open-Loop System";
          "Poles: Closed-Loop System";
          "Poles: Closed-Loop System";
          "Poles: Closed-Loop System"];

Fig_names = ["Open-loop";
             "mu = 10";
             "mu = 100";
             "mu = 10^10"];

%% Open-loop

K_final = zeros(2,3);
pos_polos(F, E_F, G, E_G, H, K_final, 1000, Fig_names(1), Titles(1));
pause

%% EigenValues (Feedback System)

for i=1:length(mu)

    % Robust LQR
    for j = 1:Horizonte
        [L(:,:,j), K(:,:,j), P(:,:,j+1)] = robust_lqr(F, G, E_F, E_G, H, Q, R, P(:, :, j), mu(i));
    end
    
    fprintf('Valores de P, K e L para mu = %f\n', mu(i))
    disp(P(:,:,Horizonte+1));
    disp(K(:,:,Horizonte));
    disp(L(:,:,Horizonte));
    fprintf('------------------------------------------\n')
    
    % Robust Control for Unceartain Systems
    [x_inc_rob, u_inc_rob] = closed_loop(K, L, xo, Horizonte, Fig_names(i+1));

    if i == length(mu)
    
        % Nominal Control for Unceartain Systemss
        [Gain, X, E] = dlqr(F, G, Q, R, zeros(3,2));

        [x_inc_nom,u_inc_nom,J_mean] = ...
        closed_loop_via_lqr_nom(F, E_F, G, E_G, H, Q, R, Gain, xo, Horizonte);

        pos_polos(F, E_F, G, E_G, H, K(:, :, Horizonte), 1000, Fig_names(i+1), Titles(i+1));

        fprintf('Fim\n')

    else

        pause
        % close all
        pos_polos(F, E_F, G, E_G, H, K(:, :, Horizonte), 1000, Fig_names(i+1), Titles(i+1));

    end
end
