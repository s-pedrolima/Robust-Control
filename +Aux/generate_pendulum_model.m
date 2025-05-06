function [F, G, H, D] = generate_pendulum_model(g, I_n, M, m, l, ...
b, Xi, k_m, Ts)
% Auxiliary variables
Delta_f = -m * l^2 * M - I_n * (m + M);
a11 = (Xi * (m + M))/Delta_f;
a12 = (-b * m * l)/Delta_f;
a13 = (-m * g * l * (m + M))/Delta_f;
a21 = (-Xi * m * l)/Delta_f;
a22 = (b * (I_n + m * l^2))/Delta_f;
a23 = (g * (m * l)^2)/Delta_f;
b1 = (k_m * m * l)/Delta_f;
b2 = (-k_m * (I_n + m * l^2))/Delta_f;
% Continous-time system matrices
A = [a11 a12 a13;
    a21 a22 a23;
    1 0 0];
B = [b1;
    b2;
    0];
n_size = size(A,1);
% Augmented system
C_a = [0 -1 0];
A_aug = [A, zeros(n_size,1);
        C_a, 0];
B_aug = [B;
        0];
n_aug_size = size(A_aug,1);
C_aug = eye(n_aug_size);
D_aug = 0;
% State space system in continuous time
pendulo_css = ss(A_aug, B_aug, C_aug, D_aug);
% State space system in discrete time
pendulo_dss = c2d(pendulo_css, Ts);
% Discrete-time system matrices
F = pendulo_dss.a;
G = pendulo_dss.b;
H = pendulo_dss.c;
D = pendulo_dss.d;
end
