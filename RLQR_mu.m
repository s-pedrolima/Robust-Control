%% RLQR Numerical Example

clc;
clearvars;
close all; 
import Controllers.*
import Aux.*

%% System Configuration

% System matrices
F = [1.1, 0, 0;
     0, 1, 0;
    -1, 1, .9];

G = [0, 1;
     1, 1;
    -1, 1];

% Include unceartainties? Y(1)/N(0)
a = 1;

% Unceartainties matrices
H = [.1;
     1.5;
    -1];

E_F = a*[.4, .5, -.6];
E_G = a*[.4, -.4];

% Initial states
xo = [1; -1; .5];

% Penalty parameters: 
mu = [1e2, 1e5, 1e10, 1e12]; % mu = 1ex > 0 

% Loop length
N = 70;

% Covariance matrices
[n_size, m_size] = size(G);
Q = eye(n_size, n_size);          % State weighting matrix
R = eye(m_size, m_size);          % Control weighting matrix
K = zeros(m_size, n_size, N);     % Gain matrix
L = zeros(n_size, n_size, N);     % L matrix
P = zeros(n_size, n_size, N+1);   % Ricatti matrix
P(:, :, 1) = eye(n_size);         % Initial P value

% Figures settings
Titles = ["Poles: Open-Loop System";
          "Poles: Closed-Loop System";
          "Poles: Closed-Loop System";
          "Poles: Closed-Loop System";
          "Poles: Closed-Loop System"];

Fig_names = ["Open-loop";
             "mu = " + mu(1);
             "mu = " + mu(2);
             "mu = " + mu(3)
             "mu = " + mu(4)];

%% Open-loop

K_final = zeros(2,3);
pos_polos(F, E_F, G, E_G, H, K_final, 1000, Fig_names(1), Titles(1));
disp("System Paused. Press any key to continue...");
pause

%% EigenValues (Feedback System)

for i=1:length(mu)

    % Robust LQR
    for j = 1:N
        [L(:,:,j), K(:,:,j), P(:,:,j+1)] = robust_lqr(F, G, E_F, E_G, H, Q, R, P(:, :, j), mu(i));
    end

    disp("------------------------------------------");
    disp("P, K and L values for mu = 10^" + log10(mu(i)));
    disp(P(:,:,N+1));
    disp(K(:,:,N));
    disp(L(:,:,N));
    disp("------------------------------------------");
    disp("Cost and EF+EG*K values for mu = 10^" + log10(mu(i)));
    disp(xo'*P(:,:,N+1)*xo);
    disp(E_F+E_G*K(:,:,N));
    disp("------------------------------------------");
   
    % Robust Control for Unceartain Systems
    [x_inc_rob, u_inc_rob] = closed_loop(K, L, xo, N, Fig_names(i+1));

    if i == length(mu)
    
        % Nominal Control for Unceartain Systemss
        [Gain, X, E] = dlqr(F, G, Q, R, zeros(3,2));

        [x_inc_nom,u_inc_nom,J_mean] = ...
        closed_loop_via_lqr_nom(F, E_F, G, E_G, H, Q, R, Gain, xo, N);

        pos_polos(F, E_F, G, E_G, H, K(:, :, N), 1000, Fig_names(i+1), Titles(i+1));

        disp("End")

    else
        disp("System Paused. Press any key to continue...");
        pause

        % Optional: close all figures and clear command window
        % close all
        % clc

        pos_polos(F, E_F, G, E_G, H, K(:, :, N), 1000, Fig_names(i+1), Titles(i+1));

    end
end
