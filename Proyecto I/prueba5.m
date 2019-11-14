function [x0, z0, ban, iter] = prueba5(A, b, c)
    
    % Inicializamos la salida
    x0 = []
    z0 = []
    ban = 0
    iter = 0
    bandera = 0;
    [m, n] = size(A)
    
    % Conjuntos de �ndice iniciales
    N = 1:n
    B = n+1:m+n
    
    % Costos b�sicos y no b�sicos
    cN = c
    cB = zeros(1, m)
    
    % Matriz b�sica, no b�sica y matriz extendida iniciales
    AN = A
    AB = eye(m)
    lambda = (AB' \ cB')'
    rN = lambda * AN - cN
    h = AB \ b'
    
    % M�todo de mayor descenso
    while bandera == 0
        if ~all(rN <= 0) == 1
                       
            %Buscamos el �ndice del que entra.
            [maxi,e] = max(rN)
            
            %Encontrar el �ndice del que sale.
            h = AB \ b'
            He = AB \ AN(:, e)
            AN(:,N == e)
            [mini,s] = min(h ./ He)
            
            % Revisamos si es no acotado.
            if all(He <= 0)
                xo = []
                zo = []
                ban = 1
                iter 
                return
            end
            N(e) = B(s)
            B(s) = e
            
            AN(:, e) = AB(:, s)
            AB(:, s) = He
            
            aux = cB(s)
            cB(s) = cN(e)
            cN(e) = aux 
            
            lambda = (AB' \ cB')'
            %lam = cB * inv(AB)
            rN = lambda * AN - cN
            
            iter = iter + 1  
        else
            xB = AB \ b'
            zo = lambda * b'
            return
        
    end
    
    
    
    
    
    
end