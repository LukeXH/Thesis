clear

% This code is an exploration of geodesics on an R2 euclidan manifold, M,
% with the euclidian metric, I_2x2

n = 2;
x = sym('x', [n,1]);
v = sym('v', [n,1]);

M = x; % Manifold parameter "vector"
g = eye(n);% Metric on M
TM = [x;v]; % TM parameter "vector"

% Now to calculate the Sasaki metric
gam1 = Christoffel_1st(g,x);
gam2 = Christoffel_2nd(g,x);
% For the sake of argument


A = zeros(n);

for j = 1:n
    for k = 1:n
        A_tmp1 = gam2(:,:,j)*diag(v);
        A_tmp2 = gam2(:,:,k)*diag(v);
        A(j,k) = 
    end
end
% B = 
g_bar = [g, 0; 0, g] + [A, B; B.', zeros(size(g))];

