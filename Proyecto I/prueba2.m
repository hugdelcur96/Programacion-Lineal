function [x0, z0, ban, iter] = prueba2(A, b, c)
    
    [m, n] = size(A);   % número de ecuaciones
    AB = eye(m);         % matriz de básicas, identidad
    AN = A;              % matriz de no basicas, A
    cB = zeros(1, m);   % costos basicas
    cN = c;             % costos no basicas
    lambda = cB * inv(AB);
    rB = zeros(1, m);
    rN = lambda * AN - cN;
    B = linspace(n+1, m+n, m);   %conjunto básicas
    N = linspace(1, n, n);   % conjunto no básicas
    
    x0 = zeros(n,1);
    z0 = cB * inv(AB) * b';
    ban = 1;
    iter = 1;
    
    fprintf("Tabla inicial")
    T = [ inv(AB) * AN, AB, inv(AB) * b'; cB * inv(AB) * AN - cN, zeros(1, m), z0]
    while T(m+1, N) > 0
        
        
        
        iter
        if T(m+1, N) <= 0
            x0 = inv(AB) * b'
            z0 = lambda * b'
            ban = 0
            break
        else
            [mx, e] = max(T(m+1, 1:n+m))
            
            if T(1:m, e) <= 0
                x0
                ban = 1
                break
            else
                [mn, sal] = min(T(1:m, n+m+1) ./ T(1:m, e))
                s = sal + n
                pivote = T(sal, e)
                T(sal, :) = T(sal, :) / pivote;
                
                for i = 1:m+1
                    if i ~= sal
                        T(i, :) = T(i, :) - T(i, e) * T(sal, :);
                    end
                end
                
            end
            
            B(B == s) = e
            N(N == e) = s
            iter = iter + 1;
            T
            
        end
    end
      
    
end






