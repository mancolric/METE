function Problema2_2026Extra()

    %Datos:
    m       = 1.2;
    c       = 0.1;
    k       = 1.5;
    
    %Condici'on inicial:
    y0      = [ 1.2; 2.1 ];
    yp0     = [ 0.3; -0.4];
    
    %Ecuaci'on a integrar:
    function ypp = fun(t, y, yp)
        
        ypp     = [ -k/m * (2*y(1)-y(2)) - c/m * (2*yp(1)-yp(2));
                    -k/m * (y(2)-y(1)) - c/m * (yp(2)-yp(1)) ];
                
    end

    %Llamamos al RK:
    [tv, ym, ypm] = RungeKuttaNystromE(@fun, 0.0, y0, yp0, 'RKN4', 0.1, 20.0);
    
    %Representamos y1 e ydot2:
    figure()
    plot(tv, ym(1,:), 'b')
    hold on
    plot(tv, ypm(2,:), 'g')
    xlabel('t')
    legend('y_1', 'dy_2/dt', 'location', 'best')
    
end
    