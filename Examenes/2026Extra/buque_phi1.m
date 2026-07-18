function [phi_k, beta_k] = buque_phi1(k, xhat) 
%Esta función es v'alida tanto si xhat es un escalar como si es un vector.

    %Calculamos beta:
    if k==1
        beta_k  = 0.0;
    else
        beta_k  = (4*k-5)*pi/4 + exp(-(4*k-5)*pi/2);
    end
    
    %Calculamos funci'on:
    phi_k       = cosh(beta_k*xhat)/cosh(beta_k) + cos(beta_k*xhat)/cos(beta_k);
    
end
