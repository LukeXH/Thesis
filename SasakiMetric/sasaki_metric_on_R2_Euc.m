clear

% This code is an exploration of geodesics on an R2 euclidan manifold, M,
% with the euclidian metric, I_2x2

x = sym('x', [2,1]);
v = sym('v', [2,1]);

M = x; % Manifold parameter "vector"
g = ones(2);% Metric on M
TM = [x;v]; % TM parameter "vector"

% Now to calculate the Sasaki metric

A = 
B = 
g_bar = [g, 0; 0, g] + [A, B; B.', zeros(size(g))];

