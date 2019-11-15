function [xo, zo, ban, iter] = prueba5(A, b, c)
    
    % Inicializamos la salida
    xo = [];
    zo = [];
    ban = 0;
    iter = 0;
    bandera = 0;
    [m, n] = size(A);
    
    % Conjuntos de �ndice iniciales
    N = 1:n;
    B = n+1:m+n;
    
    % Costos b�sicos y no b�sicos
    cN = c;
    cB = zeros(1, m);
    
    % Matriz b�sica, no b�sica y matriz extendida iniciales
    AN = A;
    AB = eye(m);
    lambda = (AB' \ cB')';
    rN = lambda * AN - cN;
    h = AB \ b';
    
    if ~all(b < 0) == 0
        xo = [];
        zo = [];
        ban = -1;
        iter;
        return
    end
    
    % M�todo de mayor descenso
    while bandera == 0
        if ~all(rN <= 0) == 1
            %Buscamos el �ndice del que entra.
            [maxi, t] = max(rN);
            e = N(t);
            %Encontrar el �ndice del que sale.
            h = AB \ b';
            He = AB \ AN(:, t);
            optTest = He > 0;       % cu�ntos valores son positivos de la columna pivote
            if sum(optTest) > 0;
                mrt = find(He > 0); % �ndices de la columna pivote que son positivos
                hs_div = h(mrt) ./ He(mrt);     % divisi�n de h/hl
                [mini, r] = min(hs_div);     % valor m�nimo e �ndice de las divisiones
                r = mrt(r);          % este ya es el �ndice de entrada final
                s = B(r);
                
                N(t) = B(r);
                B(r) = e;
                
                aux1 = AB(:,r);
                AB(:, r) = AN(:,t);
                AN(:, t) = aux1;
                
                aux2 = cB(r);
                cB(r) = cN(t);
                cN(t) = aux2; 

                lambda = (AB' \ cB')';
                %lam = cB * inv(AB)
                rN = lambda * AN - cN;

                iter = iter + 1; 
            else % cuando es no acotado
                xo = [];
                zo = [];
                ban = 1;
                iter;
                return
            end
        else
            xB = AB \ b';
            xo(B) = xB;
            xo(N) = 0;
            xo = xo(1:n)';
            zo = lambda * b';
            ban = 0;
            iter;
            return
    end
    
end