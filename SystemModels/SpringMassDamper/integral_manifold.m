function integral_manifold

m = 1;      % Mass
c = 1;      % Damping
k = 1;      % Spring Constant
w = .731/(2*pi);  % Frequency
w = 1/(2*pi);
% w = .1;

g = 0.2;    % Input gain

fcn_phys = @(q,u) [q(2);...
                   (u - c*q(2) - k*q(1))/m];
               
fcn_ctrl_kick = @(t, q) g * sign(sin(2*pi*w*t));
fcn_ctrl_push = @(t, q) g * sin(2*pi*w*t);
fcn_ctrl_opt1 = @(t, q) g * sign(q(2));
fcn_ctrl_opt2 = @(t, q) g * sign(q(2))*exp(-1000*q(2).^2);
fcn_ctrl_opt3 = @(t, q) g * sign(q(2))*(q(2) < g/c);
fcn_ctrl_circ = @(t, q) sign(q(2))*sqrt(g^2 - q(1)^2 - q(2)^2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Select Control Strat
switch 2
    case 1
       fcn_u = fcn_ctrl_kick; 
    case 2
       fcn_u = fcn_ctrl_push;
    case 3
       fcn_u = fcn_ctrl_opt1;
    case 4
       fcn_u = fcn_ctrl_opt2;
    case 5
       fcn_u = fcn_ctrl_opt3;
    case 6
        fcn_u = fcn_ctrl_circ;
    otherwise
        fcn_u = @(t, q) 0;
end


fcn_sim = @(t,x) fcn_phys(x, fcn_u(t,x));

[t,q] = ode45(fcn_sim, [0,3/w], [0;.00001]);
% Calculate the control cmds sent
u_thang = zeros(1, size(q,1));
% Recreate the U-command
for i=1:size(q)
    u_thang(i) = fcn_u(t(i),q(i,:)');
end
dt = [diff(t)', diff(t((end-1):end))];

%%% PLOT %%%
%%% NOW WITH ONE MORE DIMENSION
% Plot the vector graphgs for generated
% power, power cost and total
u_s = linspace(min(u_thang), max(u_thang), 4);
x_s = linspace(min(q(:,2)), max(q(:,2)), 4);
v_s = linspace(min(q(:,2)), max(q(:,2)), 4);
[U,X,V] = meshgrid(u_s,x_s,v_s);
P_gen_u = 0*U;
P_gen_x = U;
P_gen_v = 0*V;
P_cost_u = 0*U;
P_cost_x = -2*c*V;
P_cost_v = 0*X;
P_total_u = P_gen_u + P_cost_u;
P_total_x = P_gen_x + P_cost_x;
P_total_v = P_gen_v + P_cost_v;


figure(17301)
clf
quiver3(U,X,V, 0*U, V, (U - k*X - c*V)/m)
hold on
plot3( u_thang,q(:,1), q(:,2), 'r')
quiver3(u_thang', q(:,1), q(:,2),...
        0*u_thang', q(:,2), (u_thang' - k*q(:,1) - c*q(:,2))/m, 'Color', [.5, 0, .5])
xlabel('u')
ylabel('x')
zlabel('v')
title('Natural Dynamics')
axis square


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Looking at Lie brackets and such now
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Vector field functions
% Passive dynamics
fcn_g_1 = @(q,u) [q(2);...
                  -(c*q(2) + k*q(1))/m];
% Force input
fcn_g_2 = @(q,u) [0;...
                  u/m];
              
% Symbols
x = sym('x',[2,1]);
u = sym('u');
g_1 = fcn_g_1(x,u);
g_2 = fcn_g_2(x,u);
% Lie Bracket
Lie = @(g1, g2) [diff(g2,x(1)), diff(g2,x(2))] * g1 - ...
                [diff(g1,x(1)), diff(g1,x(2))] * g2;  

G_orig = [g_1, g_2];
G_new = Lie(g_1,g_2);
G = [G_orig, G_new];
% Generate all possible Lie Brackets
for i=2:4
    G_old = G_new;
    G_new = [];
    for n = 1:size(G_old,2)
        G_new =[G_new, Lie(g_1, G_old(:,n))];
        G_new =[G_new, Lie(g_2, G_old(:,n))];
    end
    G = [G, G_new];
end

disp(G)
disp('We just need')
disp(G(:,1:3))
a = sym('a',[3,1]);
spanG = G(:,1:3)*[1;a(2:3)];
disp('Span of G')
disp(spanG)


%%% Plots
figure(17302)
clf
[X2,V2] = meshgrid(x_s,v_s);
subplot(1,3,1)
quiver(X2,V2,V2,-k/m*X2 - c/m*V2)
xlabel('x')
ylabel('v')
title('Passive')
subplot(1,3,2)
quiver(X2,V2,zeros(size(X2,1),size(X2,2)),ones(size(X2,1),size(X2,2))/m)
xlabel('x')
ylabel('v')
title('Active (Scales by u)')
subplot(1,3,3)
quiver(X2,V2,-ones(size(X2,1),size(X2,2))/m,c/m*ones(size(X2,1),size(X2,2)))
xlabel('x')
ylabel('v')
title('Lie Bracket (Scales by u)')

%%% Integral manifold in the (t,x)-space
figure(17303)
clf
plot3(t,q(:,1),q(:,2))
xlabel('time')
ylabel('x')
zlabel('v')
title('Integral manifold in (t,x)-space')
%%% Let's do some foliation!


%%% Integral Manifold slices for different frequency sinusoidal inputs
n_data = 100;
sample_times = linspace(.1,15,50);
tspan = [0,sample_times(end)];
data_x = zeros(n_data, size(sample_times, 2));
data_v = zeros(n_data, size(sample_times, 2));
data_u = zeros(n_data, size(sample_times, 2));
for i = 1:n_data
%    sol = ode45(@(t,x) fcn_phys(x,10*(i/n_data)*g*sin( .1*2*pi*w*t)), tspan, [0;0]);
   sol = ode45(@(t,x) fcn_phys(x,g*sin(1*pi*(2*(i/n_data-.5))*t)), tspan, [0;0]);
%    foo = generate_rand_func(2*(rand(12,1)-0.5));
%    sol = ode45(@(t,x) fcn_phys(x,g*foo(t/tspan(end))), tspan, [0;0]);
   for n = 1:size(sample_times, 2)
    disp(sprintf('%.d, %.d',i,n))
    tmp = deval(sol,sample_times(n));
    data_x(i,n) = tmp(1);
    data_v(i,n) = tmp(2);
    data_u(i,n) = g*sin(1*pi*(i/n_data)*sample_times(n));
   end
end
figure(17304)
clf
for n = 1:size(sample_times, 2)
    plot3(sample_times(n)*ones(n_data,1), data_x(:,n), data_v(:,n),'.')
    hold on
end
for n = 1:n_data
   plot3(sample_times, data_x(n,:), data_v(n,:),'Color',[.8,.8,.8])%)[.5,.5*n/n_data,.5])
   hold on
end
xlabel('time')
ylabel('x')
zlabel('v')
title('Monte-Carlo Integral Manifold in (t,x)-space')
% axis equal

figure(17305)
clf
for n = 1:n_data
    plot3(data_u(n,:), data_x(n,:), data_v(n,:),'Color',[.6,.6,.6])
    hold on
end
% for n = 1:n_data
%    plot3(sample_times, data_x(n,:), data_v(n,:),'Color',[.8,.8,.8])%)[.5,.5*n/n_data,.5])
%    hold on
% end
xlabel('u')
ylabel('x')
zlabel('v')
title('Monte-Carlo Integral Manifold in (u,x)-space')
% axis equal


end