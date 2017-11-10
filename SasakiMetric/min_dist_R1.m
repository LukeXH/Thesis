clear

% This script looks at the total momentum change in R1 between two points
% in the state space.

% Parameters
q_start = [1;1];
q1 = [0;0];
q2 = [2;0];

% Relaxing function, only moves velocities up and down
da = @(v) [0; (v(1:end-2) + v(3:end))/2-v(2:end-1); 0];
% Get elements summed
f_sum = @(v,dt) dt/2*sum(v(1:end-1)+v(2:end));
% Constraint function, moves points around, points are (t, v)
ds = @(s,x_goal) [0, 0;...
                  cumsum(ones(size(s,1)-1,1)*(x_goal^2 - f_sum(s(:,2),s(2,1))^2)), [ones(size(s,1)-2,1)*(x_goal - f_sum(s(:,2),s(2,1))); 0]];... 
% Time minimization

% Finite elements
n = 10;
velo = [q_start(2), zeros(1,n-1)];
t = linspace(0,1,n);

% simulate
m = 1000;
e = zeros(m+1,1);
x = [t', velo'];
for i=1:m
    e(i) = q1(1)-q_start(1) - f_sum(x(:,2),mean(diff(x(:,1))));
    x = x + .1*ds(x,q1(1)-q_start(1)) + 0.01*[zeros(n,1), da(x(:,2))];
end
e(i+1) = q1(1)-q_start(1) - f_sum(x(:,2),mean(diff(x(:,1))));

figure(30000)
clf
subplot(2,1,1)
plot(e)
ylabel('error')
subplot(2,1,2)
plot(x(:,1),x(:,2))
ylabel('v')