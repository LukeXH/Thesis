function G = Christoffel_1st(metric, coordinates)
%CHRISTOFFEL_1ST Summary of this function goes here
%   Detailed explanation goes here
K = size(coordinates,1);
I = size(metric,1);
J = size(metric,2);
G = zeros(K,I,J);

for k = 1:K
    for i = 1:I
        for j = 1:J
            G(k,i,j) = .5*( diff(metric(k,i), coordinates(j)) + ...
                            diff(metric(k,j), coordinates(i)) - ...
                            diff(metric(i,j), coordinates(k)));
        end
    end
end
end

