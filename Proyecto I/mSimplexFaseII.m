function [x0, z0, ban, iter] = mSimplexFaseII(A, b, c)

    %clear all
    %clc

    [m, n] = size(A);
    
    if m ~= size(b)
        error('El número de renglones de A y la dimensión de b no coinciden.')
    elseif n ~= size(c)
        error('El número de columnas de A y la dimensión de c no coinciden,')
    elseif b < 0
        error('El vector b tiene que ser mayor o igual que cero.')
    end
    
    %%%%% Se agregan variables de holgura %%%%%
    
    A_simplex = zeros(m + 1, n + m + 1);
    b_simplex = zeros(m + 1, 1);    
    c_simplex = zeros(n + m + 1, 1);
    
    A_simplex(1:m, 1:n) = A;
    A_simplex(1:m, n+1:n+m) = eye(m);
    b_simplex(1:m) = b;
    c_simplex(1:n) = c;
    
    A_simplex(m + 1, 1:n+m+1) = -c_simplex;
    A_simplex(1:m, n+m+1) = b;
    
    
    % Empezamos el algoritmo del método simplex por la regla de mayor descenso.
    
    iter = 0;
    
    MaxNumIter = 100;
    
    for iteraciones = 1:MaxNumIter
        
        fin = A_simplex(m+1, 1:n+m) > 0;
        
        if fin <= 0
            break
        end
        
        [a, e] = max(A_simplex(m+1, :)); % a es el máximo, e es el índice de entrada.
        
        Xre = A_simplex(:, n+m+1) ./ A_simplex(:, e);
        i = Xre <= 0;
        d = Xre;
        d(i) = inf;
        
        [q, s] = min(d); % q es el valor mínimo, s es el índice de salida.
        
        A_simplex(s, 1:n+m+1) = A_simplex(s, :) / A_simplex(s, e);
        
        for i = 1:1:m+1
            if i ~= s
                A_simplex(i, :) = A_simplex(i, :) - A_simplex(i, e) * A_simplex(s, :);
            end
        end
        
        iter = iter + 1;
    end
    
end

