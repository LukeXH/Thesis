function g_bar = calcSasakiMetric( metric, x, v )
%CALCSASAKIMETRIC CALCULATES THE SASAKI METRIC FOR A GIVEN metric WITH THE
%POSITIONAL COORDINATES x AND VELOCITY COORDINATES v
%   Detailed explanation goes here

% Now to calculate the Sasaki metric
gam1 = christoffel1st(metric, x);
gam2 = christoffel2nd(metric, x);

%
n = size(metric, 1);

A = sym('x')*zeros(n);
B = sym('x')*zeros(n);

for j = 1:n
    for k = 1:n
        A(j,k) = v.'* (gam2(:,:,j)*metric*gam2(:,:,k)) *v;
        B(j,k) = gam1(k,:,j) * v;
    end
end

g_bar = [metric, zeros(n); zeros(n), metric] + [A, B; B.', zeros(n)];

end