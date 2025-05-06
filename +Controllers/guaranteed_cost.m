function [K, P] = guaranteed_cost(F, G, EF, EG, H, Q, R, P, epsilon)
    % Function Guaranteed Cost
    inv_aux = eye(size(H,2)) + epsilon * (H' * P * H);

    aux1 = (F' * P * F) + epsilon * ((P * H) * (inv_aux \ (H' * P)));
    
    aux2_1 = (F' * P * G) + (1/epsilon) * (EF' * EG);
    aux2_2 = (G' * P * G) + R + (1/epsilon) * (EG' * EG);
    aux2_3 = (G' * P * F) + (1/epsilon) * (EG' * EF);
    
    aux2 = aux2_1 * (aux2_2 \ aux2_3);
    
    aux3 = (1/epsilon) * (EF' * EF) + Q;
    
    P = aux1 - aux2 + aux3;
    K = - ( (G' * P * G) + R + (1/epsilon) * (EG' * EG) ) \ ...
          ( (G' * P * F) + (1/epsilon) * (EG' * EF) );

end