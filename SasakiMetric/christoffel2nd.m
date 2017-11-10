function G2 = christoffel2nd(metric, coordinates)
%CHRISTOFFEL_2ND Calculates Christoffel symbols of the second kind.
%   Detailed explanation goes here
G = christoffel1st(metric, coordinates);
g_inv = inv(metric);
n = size(g_inv, 2);

G2 = sym('x',[1,1])*zeros(size(G));
for m = 1:n
    for i=1:n
        for j=1:n
            G2(m,i,j) = g_inv(m,:)*reshape(G(:,i,j), [n,1]);
        end
    end
end
end
