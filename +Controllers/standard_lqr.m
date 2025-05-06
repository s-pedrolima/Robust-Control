function [L, K, P] = standard_lqr(F, G, Q, R, Pp)

    aux_inverse = inv(R + G' * Pp * G) * G' * Pp;

    K = -aux_inverse * F;
    L = F - G * aux_inverse * F;
    P = F' * (Pp - Pp * G * aux_inverse) * F + Q;

end 
