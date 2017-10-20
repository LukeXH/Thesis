function func_handle = generate_rand_func( weights )
%GENERATE_RAND_FUNC returns a function handle to a polynomial of degree n-1
%where n is the length of the weight vector passed in. 
%enforces the domain x = [0,1]

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

% Reshape into column vector
w = reshape(weights, [n,1]);
    
t = linspace(0,1,n+1);
t = reshape(t(2:end),[n,1]);

A = ones(n,n);

x = sym('x');
P = [];
for i = 1:n
   A(:,i) = t.^i; 
   P = [P,x^i];
end

    func_handle = matlabFunction(P*(A\w));
end

