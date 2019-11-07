function [x0, z0, ban, iter] = mSimplexFaseII(A, b, c)

    %clear all
    %clc

    [m, n] = size(A);
    
    if m ~= size(b)
        error('El número de renglones de A y la dimensión de b no coinciden.')
    elseif b < 0
        error('El vector b tiene que ser mayor o igual que cero.')
    end
    
    A_simplex = zeros(m, n + m);
    
    A_simplex(1:m, 1:n) = A;

    for i = 1:m
        A_simplex(i, n + i) = 1;
    end
    
end