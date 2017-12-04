% This script seeks to look at the evolution of a pdf of the position of a
% kinematic cart. (A vehicle that only has two wheels and moves on a
% plane).
clear

% Add 
matlab_path = pwd;
matlab_path = matlab_path(1:strfind(pwd, 'MATLAB')+5);
addpath([matlab_path,'\Handy-Tools.git\GroupTheory'])

% Assume the system evolves on SE(2), and the cart is simply trying to
% drive in a straight line.
r = .5; % wheel radius
L = 1;  % Wheel base width
D = 2; % noise variance
wheel_speed = 10;
dq = @(q, w, dw) [r*w*cos(q(3)), r*w*sin(q(3)), 0]' +...
                 sqrt(D)*[r/2*cos(q(3)), r/2*cos(q(3));...
                          r/2*sin(q(3)), r/2*sin(q(3));...
                          r/L,          -r/L]*dw;
f_rand = @(t) rand(2,1);

f_ode = @(t,q) dq(q,wheel_speed,f_rand(t));

sol = ode45(f_ode,[0,10],zeros(3,1));

figure(16700)
clf
subplot(2,1,1)
plot(sol.y(1,:),sol.y(2,:))
xlabel('x')
ylabel('y')
subplot(2,1,2)
plot(sol.x,sol.y(3,:))

% Now let's do this the Fokker-Planck way, but look at the 

