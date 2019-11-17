%%%%% EJEMPLO 1 %%%%%

fprintf("\nEJEMPLO 1")

A = [1 2 1; 2 -1 3];
b = [3 4];
c = [2 3 4];

[xo, zo, ban, iter, lamo] = mSimplexDual(A, b, c)

fprintf("Fin EJEMPLO 1\n")

%%%%% EJEMPLO 2 %%%%%

fprintf("\nEJEMPLO 2")

A = [1 2 1; 2 -1 3];
b = [3 4];
c = [2 3 4];

[xo, zo, ban, iter, lamo] = mSimplexDual(A, b, c)

fprintf("Fin EJEMPLO 2\n")
