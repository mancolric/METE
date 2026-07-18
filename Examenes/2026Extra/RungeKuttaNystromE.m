%Devuelve vector fila tv=[t_0, t_1, ..., t_n] y matriz 
% ym = [y_0, y_1, ..., y_n], donde y_i = y(t_i) es un vector columna:
% fun es una función de la forma f = fun(t, y)
function [tv, ym, ypm] = RungeKuttaNystromE(fun, t0, y0, yp0, metodo, Deltat, tf)

    %Condici'on inicial:
    tn      = t0;
    yn      = y0(:);
    ypn     = yp0(:);
    tv      = [tn];
    ym      = [yn];
    ypm     = [ypn];
    
    %Cargamos coeficientes:
    [A,Abar,b,bbar,c,s]     = CoefsRK2(metodo);
    
    %Matriz para guardar derivadas en las etapas intermedias:
    fm          = zeros(length(yn), s);
    
    %Bucle hasta llegar a tf:
    while tn<tf
        
        %Derivada en t^n:
        fm(:,1)     = fun(tn, yn, ypn);
        
        %Paso de tiempo:
        tnp1        = tn + Deltat;
        if tnp1>tf
            tnp1    = tf;
        end
        Deltatn     = tnp1-tn;
        
        %Recorremos etapas:
        for ii=2:s
            t_ii        = tn + c(ii)*Deltatn;
            S           = 0.0;
            Sbar        = 0.0;
            for jj=1:ii-1
                S      = S + A(ii,jj)*fm(:,jj);
                Sbar   = Sbar + Abar(ii,jj)*fm(:,jj);
            end
            y_ii        = yn + Deltatn*c(ii)*ypn + Deltatn^2*Sbar;
            yp_ii       = ypn + Deltatn*S;
            %y
            fm(:,ii)    = fun(t_ii, y_ii, yp_ii);
        end
                
        %Soluci'on en t^(n+1):
        S           = 0.0;
        Sbar        = 0.0;
        for jj=1:s
            S       = S + b(jj)*fm(:,jj);
            Sbar    = Sbar + bbar(jj)*fm(:,jj);
        end
        ynp1        = yn + Deltatn*ypn + Deltatn^2*Sbar;
        ypnp1       = ypn + Deltatn*S;
        
        %Guardamos soluci'on:
        tv          = cat(2, tv, tnp1);
        ym          = cat(2, ym, ynp1);
        ypm         = cat(2, ypm, ypnp1);
        
        %Actualizamos variables:
        tn          = tnp1;
        yn          = ynp1;
        ypn         = ypnp1;
        
    end
    
end
