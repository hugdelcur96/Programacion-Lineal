function [x0, z0, ban, iter] = prueba(A, b, c)
    
    [m, n] = size(A);   % número de ecuaciones
    AB = eye(m);         % matriz de básicas, identidad
    AN = A;              % matriz de no basicas, A
    cB = zeros(1, m);   % costos basicas
    cN = c;             % costos no basicas
    rN = - cN;
    rB = -cB;
    
    T = [       inv(AB) * AN,           AB,            inv(AB) * b';
          cB * inv(AB) * AN - cN,   zeros(1, m),    cB * inv(AB) * b']
    
    if factibilidad(T) == 1
        if optimalidad(T, m) == 1
            print("Tableau óptimo")
        else
            while optimalidad(T, m) ~= 1
                [valorMax1, e] = max(T(2, 1:end-1-m));      % e es el índice de salida
                colPivote = T(1:end-1, end) ./ T(1:end-1, e);
                positivo = colPivote >= 0;    % vemos qué entradas son positivas
                a = positivo .* colPivote;  % aqui guardamos solo los valores positivos y los negativos se vuelven cero
                a(a==0) = inf;
                [valorMin, s] = min(a);      % s es el índice de salida
                filaPivote = T(s, :);
                pivote = T(s, e);
                
                for i = 1:m+1           % hacemos Gauss
                    if i == s
                        T(i, :) = T(i, :) / pivote;
                    else
                        T(i, :) = T(i, :) - T(i, e) * T(s, :) / pivote;
                    end
                end
                T
            end
        end
    else
        print("Tableau no óptimo")
    end
end

function [f] = factibilidad(T)
    
    fact = T(1:end-1, end) >= 0;
    
    if fact
        f = 1;      % el tableau es factible
    else
        f = 0;      % el tableau no es factible
    end
    
end

function [o] = optimalidad(T, m)
    
    opt = T(2, 1:end-1-m) <= 0;
    
    if opt
        o = 1;      % la base es óptima
    else
        o = 0;      % la base no es óptima
    end
    
end


