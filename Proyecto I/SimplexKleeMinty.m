t = zeros(1, 8);
iteraciones = zeros(1, 8);
m = 3:10

for i = 1:8
    [A, b, c] = generalKleeMinty(i + 2);
    tic
    [xo, zo, ban, iter] = mSimplexFaseII(A, b, c);
    iteraciones(i) = iter;
    t(i) = toc;
end

fprintf("\t m \t\t\t numero de iteraciones \t\t cpu time\n")
fprintf('\t %0.4f \t\t\t %0.4f \t\t\t %.4f\n', [m; iteraciones; t])
