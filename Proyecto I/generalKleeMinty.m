function [A, b, c] = generalKleeMinty(m)
    
    aux1 = ones(m);
    aux2 = 2 * tril(aux1, -1);

    c = -ones(1, m);
    A = aux2 + diag(ones(m, 1));
    
    b = ones(1, m);
    for i = 2:m
        b(i) = 2^i - 1;
    end
    
end