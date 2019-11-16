function [xo, zo, ban, iter, lamo] = pruebaDual(A, b, c)
    
    xo = [];
    zo = 0;
    ban = 0;
    iter = 0;
    lamo = [];
    
    [m, n] = size(A);
    
    bandera = 0;
    
    N = 1:n;
    B = n+1:m+n;
    
    cN = c;
    cB = zeros(1, m);
    
    AN = -A;
    AB = eye(m);
    lambda = (AB' \ cB')';
    rN = lambda * AN - cN;
    h = -AB \ b';
    
    if all(c < 0) == 1;
        xo = [];
        zo = [];
        ban = -1;
        iter;
        lamo;
        return
    end
    
    while bandera == 0
        if ~all(h >= 0) == 1
            mrt = find(h < 0);
            [a, r] = min(h(mrt));
            
            r = mrt(r);
            s = B(r);
            
            H = AB \ AN;            
            Hs = H(r, :);
            optTest = Hs < 0;
            
            if sum(optTest) > 0;
                mrt = find(Hs<0);
                
                [a, t] = min(rN(mrt) ./ Hs(mrt));

                e = N(mrt(t));

                N(t) = B(r);
                B(r) = e;

                aux1 = AB(:,r);
                AB(:, r) = AN(:,t);
                AN(:, t) = aux1;
                
                aux2 = cB(r);
                cB(r) = cN(t);
                cN(t) = aux2;

                lambda = (AB' \ cB')';
                rN = lambda * AN - cN;
                h = -AB \ b';
                
                iter = iter + 1;
                
            else % no acotado
                xo = [];
                zo = [];
                ban = 1;
                iter;
                lamo = [];
            end
            
        else
            xo(B) = h;
            xo(N) = 0;
            xo = xo(1:n)';
            zo = -lambda * b';
            ban = 0;
            iter;
            lamo = -lambda';
            return
        end
        
    end
    
end
