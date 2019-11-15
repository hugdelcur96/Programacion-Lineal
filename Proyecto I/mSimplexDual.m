function [xo,zo,ban,iter,lamo] = mSimplexDual(A,b,c)
    
    % Simplex Dual
    
    % Inicializamos la salida
    x0 = []
    zo = []
    ban = 0
    iter = 0
    lamo = 0
    
    [m,n] = size(A)
    bandera = 0
    
    % Conjunto de índices iniciales
    
    N = 1:n
    B = n+1:m+n
    
    % Costos básicos y no básicos
    cN = c(N)
    cB = zeros(1, m)
    c = [cN cB]
    
    % Matriz básica, no básica y matriz extendida iniciales
    AN = - A(:, N)
    AB = eye(m)
    A = [AN AB]
    invAB = inv(AB)        % AB^(-1)
    
    xB = -invAB * b'        % soluciones básicas iniciales
    lambda = cB * invAB    % variables de holgura
    HRN = lambda * AN      % 
    rN = HRN - cN          % ahorro no básico
    
    while bandera == 0
        if ~all(xB >= 0) == 1       % Condición de SBF optima
            mrt = find(xB < 0)     % Buscamos los índices de las básicas menores que cero
            [a, r] = min(xB(mrt))  % Buscamos el mínimo de los xB
            rr = find(xB(mrt) == a)% Vemos si hay empates en los valores de las xB
            if length(rr) > 1
                r = mrt(rr)        % Hacemos el vector de los índices del mínimo de xB
                k = B(r)           % Buscamos los índices en las básicas
                [k, ~] = max(k)    % tomamos el índice mayor.
                r = find(B == k)   % 
            else
                r = mrt(r)
                k = B(r)
            end
            HRN = invAB(r,:) * AN
            mrt = find(HRN < 0)
            [a,t] = min(rN(mrt) ./ HRN(mrt))
            rr = find(rN(mrt) ./ HRN(mrt) == a)
            if length(rr) > 1
                l = N(mrt(rr))
                [l,~] = max(l)
            else
                l = N(mrt(t))
            end
            h_l = invAB * A(:,l)       % Columna pivote
            
            f = xB(r)               % Valor de variable basica que sale
            t = find(N == l)
            g = h_l(r)              % valor del pivote
            B(r) = l                % Actualizamos conjunto Basico    
            N(t) = k                % Actualizamos conjunto no básico
            AN = A(:, N)            % Actualizamos AN
            
            xB(r) = 0               % Actualizamos valor de la variable de salida
            h_l2 = h_l              % variable auxiliar
            h_l2(r) = -1            % redefinimos 
            xB = xB - (f / g) * h_l2     % Actualizamos la solución básica
            rN(t) = 0               % Actualizamos ahorro de la variable de entrada
            HRN(t) = 1
            rN = rN - a * HRN      % Actualizamos el ahorro no báscio %%%%%%%%%%%%%%%%%%
            invAB = inv(A(:, B))         % Actualizamos AB^(-1)    
            h_l(r) = g
            iter = iter + 1              % Hacemos una iteración más

        else
            x0(B) = xB
            x0(N) = 0
            z0 = c(B) * xB
            ban = 0
            iter
            rN
            B
            N
            return
        end % Fin if
    end % Fin while
    
    
    
    
    % Verificar que el dual es factible
    
end