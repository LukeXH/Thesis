clear

% This script looks at the total momentum change in R1 between two points
% in the state space.

% Parameters
q_start = [1;1];
q1 = [0;0];
q2 = [2;0];
q_goal = q1;

% Relaxing function, only moves velocities up and down
da = @(v) [0; (v(1:end-2) + v(3:end))/2-v(2:end-1); 0];
% Get elements summed
f_sum = @(v,dt) dt/2*sum(v(1:end-1)+v(2:end));
% Constraint function, moves points around, points are (t, v)
ds = @(s,x_goal) [0, 0;...
                  cumsum(ones(size(s,1)-1,1)*(x_goal^2 - f_sum(s(:,2),s(2,1))^2)), [ones(size(s,1)-2,1)*(x_goal - f_sum(s(:,2),s(2,1))); 0]];... 
% Time minimization

% Finite elements
k = 25;
velo = [q_start(2), zeros(1,k-2), q_goal(2)];
t = linspace(0,1,k);

% simulate
m = 30;
e_m = zeros(m+1,1);
x_m = [t', velo'];
for i=1:m
    e_m(i) = q_goal(1)-q_start(1) - f_sum(x_m(:,2),mean(diff(x_m(:,1))));
    x_m = x_m + .05*ds(x_m,q_goal(1)-q_start(1));% + 0.01*[zeros(k,1), da(x(:,2))];
end
e_m(i+1) = q_goal(1)-q_start(1) - f_sum(x_m(:,2),mean(diff(x_m(:,1))));

n = 2000;
e_n = zeros(n+1,1);
x_n = x_m;
for i=1:n
    e_n(i) = q_goal(1)-q_start(1) - f_sum(x_n(:,2),mean(diff(x_n(:,1))));
    tmp = ds(x_n,q_goal(1)-q_start(1));
    x_n = x_n + .1*[zeros(k,1), tmp(:,2)] + 0.01*[zeros(k,1), da(x_n(:,2))];
end
e_n(i+1) = q_goal(1)-q_start(1) - f_sum(x_n(:,2),mean(diff(x_n(:,1))));
% Get the displacement
p_m = zeros(k,1);
p_n = zeros(k,1);
for i=1:k
    p_m(i) = q_start(1) + trapz(x_m(1:i,1),x_m(1:i,2));
    p_n(i) = q_start(1) + trapz(x_n(1:i,1),x_n(1:i,2));
end

figure(30000)
clf
subplot(3,2,1)
plot(e_m)
ylabel('error')
xlabel('iteration')
subplot(3,2,3)
plot(x_m(:,1),x_m(:,2))
ylabel('v')
xlabel('time')
subplot(3,2,5)
plot(p_m, x_m(:,2))
xlabel('x')
ylabel('v')
axis equal
axis square
subplot(3,2,2)
plot(e_n)
xlabel('iteration')
subplot(3,2,4)
plot(x_n(:,1), x_n(:,2))
xlabel('time')
subplot(3,2,6)
plot(p_n, x_n(:,2))
xlabel('x')
ylabel('v')
axis equal
axis square

%%% In this section, we make some mountains
