function [A, Abar, b, bbar, c, s] = CoefsRK2(metodo)

    switch upper(metodo)
        
        case 'RKN4'
            A       = [ 0.0     0.0     0.0     0.0;
                        0.5     0.0     0.0     0.0;
                        0.0     0.5     0.0     0.0;
                        0.0     0.0     1.0     0.0 ];
            b       = [ 1/6;    1/3;    1/3;    1/6 ];
            Abar    = [ 0.0     0.0     0.0     0.0;
                        1/8     0.0     0.0     0.0;
                        1/8     0.0     0.0     0.0;
                        0.0     0.0     0.5     0.0 ];
            bbar    = [ 1/6;    1/6;    1/6;    0.0 ];
            c       = [ 0.0;    0.5;    0.5;    1.0 ];
            s       = 4;
            
        otherwise
            error('Metodo desconocido')
            
    end
            
end
