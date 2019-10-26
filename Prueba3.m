% PREGUNTA 2, PRIMERA PRACTICA DE LABORATORIO

f = -[18, 25, 10, 12, 15];

A = [1.2, 1.3, 0.7, 0, 0.5;
    0.7, 2.2, 1.6, 0.5, 1;
    0.9, 0.7, 1.3, 1, 0.8;
    1.4, 2.8, 0.5, 1.2, 0.6];

b = 40 * [4, 5, 3, 7];

Aeq = [];

beq = [];

lb = [0, 0, 0, 0, 0];

ub = [];

[x, z] = linprog(f', A, b, Aeq, beq, lb, ub)

clear all
clc