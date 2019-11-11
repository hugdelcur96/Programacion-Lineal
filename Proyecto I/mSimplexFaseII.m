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
        %%%%%% Sugerencia %%%%%%
        % Acabar el programa y regresar:
        % z0 = 0; 
        % x0 = null;
        % iter = -1;
        % ban = -1;
        % Porque el conjunto factible es vacío.
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
    N = 1:m; % Vector de no básicas.
    B = m:m+n; % Vector de básicas.
    
    
    % Empezamos el algoritmo del método simplex por la regla de mayor
    % descenso
  
    iter = 0;
    
    MaxNumIter = 100;
    
    for iteraciones = 1:MaxNumIter
        
        fin = A_simplex(m+1, 1:n+m) > 0; % Vector lógico que dice si las 
        %entradas del renglón de ahorro son menores iguales a 0 (0), mayores a 0 (1).
        c = zeros(n+m);
        if fin == c %Si todo el vector fin es cero i.e. todas las entradas
            %del vector de ahorro son menores iguales a cero, acaba el 
            % programa.
            break
        end
        
        [mx, e] = max(A_simplex(m+1, 1:n+m)); % mx es el máximo, e es el 
        %índice del máximo. Sin tener en cuenta el ahorro r0.
        
        %Definimos a h y a Hse.
        h = A_simplex(1:m, n+m+1); %Lado derecho del método simplex.
        Hse = A_simplex(1:m, e); % Columna del vector de entrada.
        
        %Creamos al vector Hse nuevo para ver qué variable sale.
        HseN = zeros(m);

        for i = 1:m %las columnas 
            if Hse(i) > 0
            HseN(i) = h(i) / Hse(i);
            end
        end
        [mn, s] = min(HseN); % mn es el mínimo, s es el índice del mínimo,
        %pero usando de referencia el vector B.
        varEntrada = e; %índice de la variable de entrada.
        varSalida = B(s); %índice de la variable de salida.
        B(s) = varEntrada; %Nueva variable básica.
        N(e) = varSalida; %Nueva variable no básica.
        
        % Dividimos el renglón de la nueva variable básica entre la 
        % entrada (varEntrada, varSalida) para volver 1.
        A_simplex(varSalida, 1:n+m+1) = A_simplex(varSalida, :) / A_simplex(varSalida, varEntrada);
        
        A_simplex(m + 1,:) = A_simplex(m + 1,:) - A_simplex(s, :) * mx;
        %Hacemos cero la entrada de la nueva variable básica.
        
        for i = 1:1:m+1
            if i ~= s
                A_simplex(i, :) = A_simplex(i, :) - A_simplex(i, e) * A_simplex(s, :);
            end
        end
        
        iter = iter + 1;
    end
    
    for i = 1:size(c, 2)
        d = logical(A_simplex(:, i));
        x0(i, 1) = A_simplex(d, end);
    end
    
end







