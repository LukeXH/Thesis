% 2 Joint Throwing Arm

matlab_path = pwd;
matlab_path = matlab_path(1:strfind(pwd, 'MATLAB')+5);
addpath([matlab_path,'\Handy-Tools.git\GroupTheory'])
addpath([matlab_path,'\Handy-Tools.git\SFDyn'])

% Joint angles
q = sym('th', [2,1]);
% Link lengths
l1 = .1;
l2 = .1;

% Set up coordinates
q1 = SE2([0,0,q(1)])*SE2([l1, 0, 0]);
q1.g = simplify(q1.g);
q2 = q1*SE2([0, 0, q(2)])*SE2([l2, 0, 0]);
q2.g = simplify(q2.g);

% Get Jacobian of end point, xy, also the precision of those variables
J = [diff(q2.xy, q(1)), diff(q2.xy, q(2))];
fun_J = matlabFunction(J);
% Sensitivity
S = [J(1,1).^2 + J(1,2).^2; J(2,1).^2 + J(2,2).^2];
fun_S = matlabFunction(S);
grad_S = simplify( [diff(S, q(1)), diff(S, q(2))] );
fun_grad_S = matlabFunction(grad_S);


% Set up grid
n = 21;
[T1, T2] = meshgrid(linspace(0, pi, n), linspace(-pi, pi, n));
ind = (-2*T1 <= T2) & ((2*pi-2*T1) >= T2);
mask = nan(n);
mask(ind) = 0;

% Set up vector fields
sol = arrayfun(fun_J, T1, T2, 'UniformOutput', false);
dxdth1 = cellfun(@(M) M(1,1), sol);
dxdth2 = cellfun(@(M) M(1,2), sol);
dydth1 = cellfun(@(M) M(2,1), sol);
dydth2 = cellfun(@(M) M(2,2), sol);
% Sensitivities (variances?)
s_dx = dxdth1.^2 + dxdth2.^2;
s_dy = dydth1.^2 + dydth2.^2;
% Sensitivity gradients
% s_dx_dth1
sol = arrayfun(fun_grad_S, T1, T2, 'UniformOutput', false);
S_dxdth1 = cellfun(@(M) M(1,1), sol);
S_dxdth2 = cellfun(@(M) M(1,2), sol);
S_dydth1 = cellfun(@(M) M(2,1), sol);
S_dydth2 = cellfun(@(M) M(2,2), sol);


%%% Plot vector field
figure(200)
clf
subplot(2,3,1)
quiver(T1(ind),T2(ind),dxdth1(ind),dxdth2(ind))
title('Vectors in X')
xlabel('th1')
ylabel('th2')
subplot(2,3,2)

surf(T1,T2,s_dx + mask)
xlabel('th1')
ylabel('th2')
subplot(2,3,3)
quiver(T1(ind),T2(ind),S_dxdth1(ind),S_dxdth2(ind))
title('Vectors in grad S')
xlabel('th1')
ylabel('th2')
subplot(2,3,4)
quiver(T1(ind),T2(ind),dydth1(ind),dydth2(ind))
title('Vectors in Y')
xlabel('th1')
ylabel('th2')
subplot(2,3,5)
surf(T1,T2,s_dy + mask)
xlabel('th1')
ylabel('th2')
subplot(2,3,6)
quiver(T1(ind),T2(ind),S_dydth1(ind),S_dydth2(ind))
title('Vectors in grad S')
xlabel('th1')
ylabel('th2')

figure(201)
clf
surf(T1,T2,s_dx + s_dy + mask)
xlabel('th1')
ylabel('th2')
title('Combined Sensitivity (Probably correct cuz you add variances)')

% figure(300)
% clf
% surf(T1, T2, 1*ind)


