function [xo, zo, ban, iter] = mSimplexFaseII(A, b, c)

    xo = [];
    zo = [];
    ban = 0;
    iter = 0;
    bandera = 0;
    [m, n] = size(A);
    
    N = 1:n;
    B = n+1:m+n;
     
    cN = c;
    cB = zeros(1, m);
    
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
    
    while bandera == 0
        if ~all(rN <= 0) == 1            
            [maxi, t] = max(rN);
            e = N(t);            
            h = AB \ b';
            He = AB \ AN(:, t);
            optTest = He > 0;
            if sum(optTest) > 0;
                mrt = find(He > 0);
                hs_div = h(mrt) ./ He(mrt);
                [mini, r] = min(hs_div);
                r = mrt(r);
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
                rN = lambda * AN - cN;
                
                iter = iter + 1; 
            else
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
