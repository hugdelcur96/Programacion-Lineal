% Una última prueba

function [x0, z0, ban, iter] = prueba4(A, b, c)
    
    % Inicializamos la salida
    x0 = [];
    z0 = [];
    ban = 0;
    iter = 0;
    
    bandera = 0;
    
    [m, n] = size(A);
    
    % Conjuntos de índice iniciales
    N = 1:n;
    B = n+1:m+n;
    
    % Costos básicos y no básicos
    cN = c(N);
    cB = zeros(1, m);
    c = [cN cB];
    
    % Matriz básica, no básica y matriz extendida iniciales
    AN = A(:, N);
    AB = eye(m);
    A = [AN AB];
    invAB = inv(AB);         % AB^(-1)
    
    xB = invAB * b';         % soluciones básicas iniciales
    lambda = cB * invAB;     % variables de holgura
    HRN = lambda * AN;       % 
    rN = HRN - cN;           % ahorro no básico
    
    
    % Método de mayor descenso
    while bandera == 0
        if ~all(rN <= 0) == 1       % Condición de SBF optima
            [~, t] = max(rN);        % índice del mayor ahorro
            l = N(t);                % ?indice en las no basicas?
            h_l = invAB * A(:, l);   % columna pivote
            optTest = h_l > 0;       % cuántos valores son positivos de la columna pivote
            if sum(optTest) > 0
                mrt = find(h_l > 0); % índices de la columna pivote que son positivos
                Xb_hl_div = xB(mrt) ./ h_l(mrt);     % división de h/hl
                [p, r] = min(Xb_hl_div);     % valor mínimo e índice de las divisiones
                rr = find(Xb_hl_div == p);   
                
                if length(rr) > 1       % el desempate cuando hay valores iguales en la división
                    r = mrt(rr);         % y se elige la variable con menor índice.
                    k = B(r);
                    [k, ~] = max(k);
                    r = find(B == k);
                else                    % Caso en donde no hay empate
                    r = mrt(r);          % este ya es el índice de entrada final
                    k = B(r);            % este ya es el índice de salida final
                end
                
                a = rN(t);               % Valor de ahorro de entrada
                f = xB(r);               % Valor de variable basica que sale
                B(r) = l;                % Actualizamos conjunto Basico
                g = h_l(r);              % valor del pivote
                N(t) = k;                % Actualizamos conjunto no básico
                AN = A(:, N);            % Actualizamos AN         
                xB(r) = 0;               % Actualizamos valor de la variable de salida
                h_l2 = h_l;              % variable auxiliar
                h_l2(r) = -1;            % redefinimos 
                xB = xB - (f / g) * h_l2;     % Actualizamos la solución básica
                rN(t) = 0;               % Actualizamos ahorro de la variable de entrada
                HRN = invAB(r, :) * AN;  
                rN = rN - (a / g) * HRN;      % Actualizamos el ahorro no báscio
                invAB = inv(A(:, B));         % Actualizamos AB^(-1)
                h_l(r) = g;                   
                iter = iter + 1;              % Hacemos una iteración más
            else        % ESTE ES EL CASO EN QUE SEA NO ACOTADO
                x0 = {};
                z0 = {};
                ban = 1;
                iter;
                return
            end 
        else
            x0(B) = xB;
            x0(N) = 0;
            z0 = c(B) * xB;
            ban = 0;
            iter;
            return
        end % Fin if
    end % Fin while
    
end % Fin funcion