function [x0, z0, ban, iter] = prueba3(A, b, c)
    
    % Inicializamos la salida
    x0 = [];
    z0 = [];
    ban = 0;
    iter = 0;
    
    [m, n] = size(A);
    
    % Conjuntos de índice iniciales
    N = 1:n;
    B = n+1:m+n;
    
    AN = A(:, N);
    AB = eye(m);
    invAB = inv(AB);
    cN = c(N);
    cB = zeros(1, m);
    
    x0(B) = invAB * b';
    w = cB * invAB;
    rN = w * AN - cN;
    
    
    
end















