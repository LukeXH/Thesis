function main

m = 1;      % Mass
c = 1;      % Damping
k = 1;      % Spring Constant
w = .731/(2*pi);  % Frequency
w = 1/(2*pi);

g = 0.02;    % Input gain

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



%%% Plots %%%
figure(12300)
clf
plot(q(:,1), q(:,2))
hold on
quiver(q(:,1),q(:,2), q(:,2), -k/m*q(:,1) - c/m*q(:,2))
quiver(q(:,1),q(:,2), q(:,2), u_thang'/m-k/m*q(:,1) - c/m*q(:,2))
title('Phase plot')
xlabel('x')
ylabel('dx')
legend('','passive','active')
axis equal

figure(12301)
clf
n_figs = 6;
% Position
subplot(n_figs,1,1)
plot(t,q(:,1))
ylabel('x')
% Velocity
subplot(n_figs,1,2)
plot(t,q(:,2))
ylabel('dx')
% Input cmd
subplot(n_figs,1,3)
plot(t, u_thang, 'r')
ylabel('u')
% Energy in system
subplot(n_figs,1,4)
plot(t, .5*k*q(:,1).^2 + .5*m*q(:,2).^2, 'm')
ylabel('E_{actual}')
% Energy generated
subplot(n_figs,1,5)
plot(t, cumsum(q(:,2)'.*u_thang.*dt), 'm');
ylabel('E_{gen}')
% Energy lost
subplot(n_figs,1,6)
plot(t, cumsum(c*q(:,2)'.^2.*dt), 'm');
ylabel('E_{cost}')
xlabel('time (sec)')
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Power integrated
% subplot(n_figs, 1, 7)
% plot(t, cumsum(q(:,2)'.*u_thang.*dt)...
%        -cumsum(c*q(:,2)'.^2.*dt), 'm')
% ylabel('E_{diff}')
% xlabel('x')
% % Energy recalc'd
% subplot(n_figs, 1, 8)
% plot(t(1:end-1), cumsum(u_thang(1:end-1) .* diff(q(:,1)')), 'k')
% ylabel('u*dx')


figure(12302)
plot(u_thang,q(:,1))
xlabel('u')
ylabel('x')
title('Velocity v. input cmd')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the vector graphgs for generated
% power, power cost and total
u_s = linspace(min(u_thang), max(u_thang), 11);
x_s = linspace(min(q(:,2)), max(q(:,2)), 11);
[U,X] = meshgrid(u_s,x_s);
P_gen_u = 0*U;
P_gen_x = U;
P_cost_u = 0*U;
P_cost_x = -2*c*X;
P_total_u = P_gen_u + P_cost_u;
P_total_x = P_gen_x + P_cost_x;

figure(12303)
clf
subplot(1,3,1)
quiver(U,X,P_gen_u, P_gen_x)
xlabel('u')
ylabel('x')
title('Power Generated')
subplot(1,3,2)
quiver(U,X,P_cost_u, P_cost_x)
xlabel('u')
ylabel('x')
title('Power Cost')
subplot(1,3,3)
quiver(U,X,P_total_u, P_total_x)
xlabel('u')
ylabel('x')
title('Power total')

%%%NOW WITH ONE MORE DIMENSION
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

figure(12304)
clf
subplot(1,3,1)
quiver3(U,X,V, P_gen_u, P_gen_x, P_gen_v)
xlabel('u')
ylabel('x')
zlabel('v')
title('Power Generated')
subplot(1,3,2)
quiver3(U,X,V, P_cost_u, P_cost_x, P_cost_v)
xlabel('u')
ylabel('x')
zlabel('v')
title('Power Cost')
subplot(1,3,3)
quiver3(U,X,V, P_total_u, P_total_x, P_total_v)
hold on
plot3( u_thang,q(:,1), q(:,2), 'r')
xlabel('u')
ylabel('x')
zlabel('v')
title('Power total')

figure(12305)
clf
quiver3(U,X,V, 0*U, V, (U - k*X - c*V)/m)
hold on
plot3( u_thang,q(:,1), q(:,2), 'r')
quiver3(u_thang', q(:,1), q(:,2),...
        0*u_thang', q(:,2), (u_thang' - k*q(:,1) - c*q(:,2))/m)
xlabel('u')
ylabel('x')
zlabel('v')
title('Natural Dynamics')



end