% Este es el mejor, solo hay que mejorar un poco la notacion y la
% explicacion.

function [x0, z0, ban, iter] = prueba3(A, b, c)
    
    % Inicializamos la salida
    x0 = []
    z0 = []
    ban = 0
    iter = 0
    
    bandera = 0
    
    [m, n] = size(A)
    
    % Conjuntos de índice iniciales
    N = 1:n
    B = n+1:m+n
    
    % Costos básicos y no básicos
    cN = c(N)
    cB = zeros(1, m)
    c = [cN cB]
    
    % Matriz básica y no básica
    AN = A(:, N)
    AB = eye(m)
    A = [AN AB]
    invAB = inv(AB)
    
    xB = invAB * b'         % soluciones básicas
    lambda = cB * invAB     % variables de holgura
    HRN = lambda * AN       
    sN = cN - HRN
    
    while bandera == 0
        if ~all(sN >= 0) == 1
            [~, t] = min(sN)
            l = N(t)
            h_l = invAB * A(:, l)
            optTest = h_l > 0
            if sum(optTest) > 0
                mrt = find(h_l > 0)
                Xb_hl_div = xB(mrt) ./ h_l(mrt)
                [p, r] = min(Xb_hl_div)
                rr = find(Xb_hl_div == p)
                
                if length(rr) > 1
                    r = mrt(rr)
                    k = B(r)
                    [k, ~] = max(k)
                    r = find(B == k)
                else
                    r = mrt(r)
                    k = B(r)
                end
                
                a = sN(t)
                f = xB(r)
                B(r) = l
                g = h_l(r)
                N(t) = k
                AN = A(:, N)
                iter = iter + 1
                xB(r) = 0
                h_l2 = h_l
                h_l2(r) = -1
                xB = xB - (f / g) * h_l2
                sN(t) = 0
                HRN = invAB(r, :) * AN
                sN = sN - (a / g) * HRN
                invAB = inv(A(:, B))
                h_l(r) = g
            else
                x0 = {}
                z0 = {}
                ban = 1
                iter
                return
            end 
        else
            x0 = xB
            z0 = c(B) * xB
            ban = 0
            iter
            return
        end
    end
end
