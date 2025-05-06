function [L, K, P] = robust_lqr(F, G, N_Fk, N_Gk, M_k, Q, R, P, mu, alpha)
    
    % Get the matrix sizes
    [n, m] = size(G); % n: number of states, m: number of control inputs
    l = size(N_Gk, 1); % l: number of rows in N_Gk (dimension of uncertainty)

    % Extended matrices
    I_r = [eye(n); zeros(l, n)]; % (n+l) x n
    G_r = [G; N_Gk];             % (n+l) x m
    F_r = [F; N_Fk];             % (n+l) x n

    % Declaring optional parameters
    if nargin < 9
        mu = 1e12; 
    end
    if nargin < 10
        alpha = 0.01/norm(mu * (M_k' * M_k), 2);
    end

    % Robustness weighting: Σ_{μ,k}
    lambda_mu = (1 + alpha) * norm(mu * (M_k' * M_k), 2); % ||.||_2
    Phi = (1/mu) * eye(n) - (1/lambda_mu) * (M_k * M_k'); % Φ(μ,λ): n x n
    Sigma_muk = blkdiag(Phi, (1/lambda_mu) * eye(l));     % Σ_{μ,k}: (n+l) x (n+l)

    % Block matrix L
    L_block = [
        zeros(n, n), zeros(n, m), zeros(n, n);
        zeros(m, n), zeros(m, m), zeros(m, n);
        zeros(n, n), zeros(n, m), -eye(n);
        zeros(n + l, n), zeros(n + l, m), F_r;
        eye(n), zeros(n, m), zeros(n, n);
        zeros(m, n), eye(m), zeros(m, n)
    ];

    % Block matrix M
    M_block = [
        P \ eye(n), zeros(n, m), zeros(n, n), zeros(n, l + n), eye(n), zeros(n, m);
        zeros(m, n), R \ eye(m), zeros(m, n), zeros(m, l + n), zeros(m, n), eye(m);
        zeros(n, n), zeros(n, m), Q \ eye(n), zeros(n, l + n), zeros(n, n), zeros(n, m);
        zeros(l + n, n), zeros(l + n, m), zeros(l + n, n), Sigma_muk, I_r, -G_r;
        eye(n), zeros(n, m), zeros(n, n), I_r', zeros(n, n), zeros(n, m);
        zeros(m, n), eye(m), zeros(m, n), -G_r', zeros(m, n), zeros(m, m)
    ];

    % Block matrix R
    R_block = [
        zeros(n, n);
        zeros(m, n);
        -eye(n);
        F_r;
        zeros(n, n);
        zeros(m, n)
    ];

    % Solve for L, K, P
    LKP = (L_block') * (M_block \ R_block);

    L = LKP(1:n, :);
    K = LKP(n+1:n+m, :);
    P = LKP(n+m+1:end, :);
end
