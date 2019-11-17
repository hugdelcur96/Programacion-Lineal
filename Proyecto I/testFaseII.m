%%%%%% EJEMPLO 1: Soluci�n �ptima %%%%%

fprintf("\nEJEMPLO 1: Soluci�n �ptima")

A = [1 0; 0 2; 3 2];
b = [4 12 18];
c = [-3 -5];
[xo, zo, ban, iter] = mSimplexFaseII(A, b, c)

fprintf("Fin EJEMPLO 1\n")


%%%%% EJEMPLO 2: Soluci�n �ptima %%%%%

fprintf("\nEJEMPLO 2: Soluci�n �ptima")

A = [1 1 2; 1 1 -1; -1 1 1];
b = [9 2 4];
c = [1 1 -4];
[xo, zo, ban, iter] = mSimplexFaseII(A, b, c)

fprintf("Fin EJEMPLO 2\n")

%%%%% EJEMPLO 3: No acotado %%%%%

fprintf("\nEJEMPLO 3: No acotado")

A = [1 -1; -1 1];
b = [1 2];
c = [-1 0];
[xo, zo, ban, iter] = mSimplexFaseII(A, b, c)

fprintf("Fin EJEMPLO 3\n")

%%%%% EJEMPLO 4: Conjunto factible vac�o %%%%%

fprintf("\nEJEMPLO 4: Conjunto factible vac�o")

A = [1 0; 0 2; 3 2];
b = [-4 -12 -18];
c = [-3 -5];
[xo, zo, ban, iter] = mSimplexFaseII(A, b, c)

fprintf("Fin EJEMPLO 4\n")
