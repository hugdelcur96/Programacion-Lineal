function [x0, z0, ban, iter] = prueba(A, b, c)
    
    m = size(A, 1);     % número de ecuaciones
    B = eye(m);         % matriz de básicas, identidad
    N = A;              % matriz de no basicas, A
    CB = zeros(1, m);   % costos basicas
    CN = c;             % costos no basicas
    
    
    
end
