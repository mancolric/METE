function gv = buque_residuo1(xhatv, A, k1, k2, eta0, av)

    %Hallamos N:
    N           = length(xhatv);
                        
    %Vector residuo precondicionado:
    gv          = zeros(N, 1);
    for ii=1:N
        
        %Valor de xhat:
        xhat_ii     = xhatv(ii);
        
        %Deflexi'on y derivada cuarta respecto de xhat en x_i:
        u_i         = 0.0;
        D4u_i       = 0.0;
        for kk=1:N
            [phik_i, betak]     = buque_phi1(kk, xhat_ii); %phi_k(x_i), beta_k
            u_i                 = u_i + av(kk) * phik_i;
            D4u_i               = D4u_i + av(kk) * betak^4 * phik_i;
        end
        
        %Altura de la ola en x_i:
        eta_i       = - eta0 * cos(pi*xhat_ii);
        
        %Residuo precondicionado:
        [phii_i, betai]     = buque_phi1(ii, xhat_ii); %phi_i(x_i), beta_i
        f_i                 = A*D4u_i + k1*(u_i - eta_i) - k2*(u_i - eta_i)^2;
        J_ii                = A*betai^4 * phii_i + ...
                                k1*phii_i - ...
                                k2*2*(u_i - eta_i)*phii_i;
        gv(ii)              = f_i/J_ii;
        
    end
    
end
