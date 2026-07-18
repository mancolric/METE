function Problema2()

    %Datos:
    m       = 1.2;
    c       = 0.1;
    k       = 1.5;
    
    %Condici'on inicial:
    y0      = [ 1.2; 2.1 ];
    yp0     = [ 0.3; -0.4];
    
    %Ecuaci'on a integrar:
    function ypp = fun(t, y, yp)
        
        
                
    end

    %Llamamos al RK:
    [tv, ym, ypm] = RungeKutta2E(   , 0.0, y0, yp0, 'RKN4',    ,    );
    
    %Representamos y1 e ydot2:
    figure()
    plot(     ,    , 'b')
    hold on
    plot(     ,    , 'g')
    xlabel('t')
    legend('y_1', 'dy_2/dt', 'location', 'best')
    
end
    
