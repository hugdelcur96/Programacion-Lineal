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
    
    
    % Empezamos el algoritmo del método simplex por la regla de Bland
    
    iter = 0;
    
    while A_simplex(m + 1, 1:n+m) > 0
        
        if A_simplex(m + 1, 1:n+m) <= 0
            break
        end
        
        indice_entrada = min(find(A_simplex(m+1, 1:n+m));
        
        if A_simplex(1:m, find(A_simplex(1:m, indice_entrada))) > 0
            indice_salida = min(A_simplex(1:m, n+m+1) ./ A(1:m, indice_entrada));
        end
        
        iter = iter + 1;
    end
    
end

%Haciendo una prueba












