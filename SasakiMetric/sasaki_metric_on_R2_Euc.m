clear

% This code is an exploration of geodesics on an R2 euclidan manifold, M,
% with the euclidian metric, I_2x2

n = 2;
x = sym('x', [n,1]);
v = sym('v', [n,1]);

M = x; % Manifold parameter "vector"
g = diag([1,x(1).^2]);%eye(n);% Metric on M
TM = [x;v]; % TM parameter "vector"

% Now to calculate the Sasaki metric
gam1 = Christoffel_1st(g,x);
gam2 = Christoffel_2nd(g,x);
% For the sake of argument


A = sym('x')*zeros(n);
B = sym('x')*zeros(n);

for j = 1:n
    for k = 1:n
        A(j,k) = v.'* (gam2(:,:,j)*g*gam2(:,:,k)) *v;
        B(j,k) = gam1(k,:,j)*v;
    end
end

g_bar = [g, zeros(size(g)); zeros(size(g)), g] + [A, B; B.', zeros(size(g))];
disp('The Sasaki Metric')
disp(g_bar)


%%% Plot metric ellipses
figure(11100)
clf
[xs,ys] = meshgrid(linspace(-1,1,5), linspace(-1,1,5));

