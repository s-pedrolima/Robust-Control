function [Lambda, P] = h_inf_ga(F, G, D, Q, Pp, gamma)
    I = eye(size(Pp,2));
    disturbance_term = (1/gamma^2) * (D * D');
    control_term = G * G';
    Lambda = I + (control_term - disturbance_term) * Pp;
    P = Q + F' * Pp / Lambda * F; 
    
    % Stability check
    Xi = gamma^2 * I - D' * P * D;
    Xi_eig = eig(Xi);
    
    if any(Xi_eig <= 0)
        error('Gamma too small: condition violated');
    end
end