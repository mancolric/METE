function Problema1_2026Extra()

    %Datos:
    N           = 6;
    A           = 2e5;
    K1          = 5e5;    
    K2          = 5e3;  
    L           = 300;
    eta0        = 7;
    h0          = 8;
    
    %Vector con xhat's:
    xhatv       = linspace(0.0, 1.0, N).';
    
    %Aproximaci'on inicial:
    av0         = zeros(N, 1);
    
    %Llamada a Anderson:
    gfun                = @(av) buque_residuo1(xhatv, A, K1, K2, eta0, av);
    [av, nIters, flag]  = Anderson1(gfun, av0, 1e-8, 100, 4);
    if flag<0
        warning('Anderson no ha convergido')
    end
    
    %Representamos ola y l'inea de base:
    xhat_plot   = linspace(-1.0, 1.0, 1000);
    x_plot      = L/2 + L/2*xhat_plot;
    eta_plot    = -eta0 * cos(pi*xhat_plot);
    u_plot      = 0*xhat_plot;
    for kk=1:N
        [phik_plot,~]   = buque_phi1(kk, xhat_plot);
        u_plot          = u_plot + phik_plot*av(kk);
    end
    
    figure()
    plot(x_plot, eta_plot, 'k')
    hold on
    plot(x_plot, -h0 + u_plot, 'b')
    xlabel('x [m]')
    grid('on')

end
