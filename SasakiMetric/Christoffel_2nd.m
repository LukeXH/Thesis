function G2 = Christoffel_2nd(metric, coordinates)
%CHRISTOFFEL_1ST Summary of this function goes here
%   Detailed explanation goes here
G = Christoffel_1st(metric, coordinates);
g_inv = inv(metric);
K = size(g_inv, 2);

G2 = zeros(size(G));
for k = 1:K
    G2(:,:,:) = G2 + g_inv(:,k)
end
end
