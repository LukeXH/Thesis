function func_handle = generate_weighted_legendre_poly( weights )
%GENERATE_LEGENDRE_POLY Returns a function handle for a legendre polynomial
%of degree n that goes from x = [-1, 1] and y = [-1, 1]
%   Detailed explanation goes here

[n1, n2] = size(weights);
if 1 == n2
    n = n1;
elseif 1 == n1
    n = n2;
else
    func_handle = @(x) 0;
    disp('Weights were not in usuable format')
    return
end

% Reshape into row vector
w = reshape(weights, [1,n]);
% w = w/sum(w);
    

x = sym('x');
P = [1; x];
for i = 1:(n-2)
    disp(i)
    P(i+2) = ( (2*i+1)*x*P(i+1) - i*P(i) )/(i+1);
end
    func_handle = matlabFunction(w*P);
end

