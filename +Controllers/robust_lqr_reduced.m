function [L,K,P] = robust_lqr_reduced(F,G,E_F,E_G,Q,R,Pp)
    % Get the matrices size
    [n, ~] = size(G);
    l = size(E_G,1);

    A_cal = [eye(n); zeros(l,n)];
    G_cal = [G; E_G];
    F_cal = [F; E_F];

    aux = (A_cal / Pp * A_cal' + G_cal / R * G_cal') \ F_cal;
    K = -(R \ G_cal') * aux;
    L = (Pp \ A_cal') * aux;
    P = F_cal' * aux + Q;
end
