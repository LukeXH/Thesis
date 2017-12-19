% Looks at 
clear


% Gaussian we're trying to approximate
u_f = 0;
s_f = 1;
f = @(x) 1/sqrt(2*pi*s_f^2) .* exp(-(x-u_f).^2./(2*s_f.^2));

% Basis guassians
p = @(x, x_i, s, w) w./sqrt(2*pi*s.^2) .* exp(-(x-x_i).^2./(2*s.^2));



% Normal Randomly sampled points
N = 20;
xs = sort(normrnd(u_f, s_f, N, 1));
[tmp_xs2, tmp_xs1] = meshgrid(xs,xs);
sig = sym('s');
weight = sym('w');
p_mat = simplify(arrayfun(@(x1,x2) p(x1,x2,sig,weight), tmp_xs1, tmp_xs2));
L_sym = ones(1,N) * (f(xs) - p_mat*(ones(N,1)/N)).^2;
L_fcn = matlabFunction(L_sym, 'Vars', {[sig; weight]});
% L_fcn = matlabFunction(L_sym);

s_w = fmincon(L_fcn, [1/N;1], -eye(2), zeros(2,1))

%%% Plotting time
y = linspace(-3, 3, 31);
figure(21300)
clf
plot(xs, zeros(N,1), 'x')
hold on
plot(y, f(y))
p_disc = zeros(N,1);
for i=1:N
   p_tmp = p(y,xs(i),s_w(1),s_w(2))/N;
   p_disc= p_disc + p_tmp;
   plot(y, p_tmp, 'g')
end
plot(y, p_disc, 'k')
