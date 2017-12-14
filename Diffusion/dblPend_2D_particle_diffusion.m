%
clear

matlab_path = pwd;
matlab_path = matlab_path(1:strfind(pwd, 'MATLAB')+5);
addpath([matlab_path,'\Handy-Tools.git\GroupTheory'])
addpath([matlab_path,'\Handy-Tools.git\SFDyn'])
addpath([matlab_path,'\Handy-Tools.git'])

% Generate equations of motion for double pendulum
% Angle of joints measured relative to the ground
dtool = DynTool;
[phi1,dphi1] = dtool.addCoord('p1');
[phi2,dphi2] = dtool.addCoord('p2');

m1 = 1;
m2 = 1;
l1 = 1;
l2 = 1;
g = 9.81;

I1 = m1 * l1^2/4;
I2 = m2 * l2^2/4;

% Neutral position is flat along the ground, position is middle of the link
pos1 = l1/2 * [ cos(phi1), sin(phi1) ];
pos2 = 2*pos1 + l2/2 * [ cos(phi2), sin(phi2) ];

vel1 = jacobian(pos1, [phi1; phi2]) * [dphi1; dphi2];
vel2 = jacobian(pos2, [phi1; phi2]) * [dphi1; dphi2];

dtool.addKE(simplify(m1/2 * dot(vel1, vel1) + I1/2 * dphi1^2));
dtool.addKE(simplify(m2/2 * dot(vel2, vel2) + I2/2 * dphi2^2));

% dtool.addPE(simplify(g * m1 * pos1(2)));
% dtool.addPE(simplify(g * m2 * pos2(2)));

sfdyn = dtool.genSFDyn;
accel_fcn = sfdyn.gen_accel_fcn;

ode_fcn = @(t, x) [x(3:4); accel_fcn(x(1:2), x(3:4), [0; 0])];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now let's try a particle swarm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
multi_ode_fcn = @(t, x, u) [x(21:40);...
                            accel_fcn(x(1:2), x(21:22), u(1:2));...
                            accel_fcn(x(3:4), x(23:24), u(3:4));...
                            accel_fcn(x(5:6), x(25:26), u(5:6));...
                            accel_fcn(x(7:8), x(27:28), u(7:8));...
                            accel_fcn(x(9:10), x(29:30), u(9:10));...
                            accel_fcn(x(11:12), x(31:32), u(11:12));...
                            accel_fcn(x(13:14), x(33:34), u(13:14));...
                            accel_fcn(x(15:16), x(35:36), u(15:16));...
                            accel_fcn(x(17:18), x(37:38), u(17:18));...
                            accel_fcn(x(19:20), x(39:40), u(19:20))];
ts = [0];
dt = .01;
x_last = zeros(40,1);
p1s = [x_last(1:2:19)];
p2s = [x_last(2:2:20)];
for i=1:1000
    rng = rand(20,1)-.5;
    tmp_sol = ode45(@(t,X) multi_ode_fcn(t,X,rng), [0, dt], x_last);
    x_last = tmp_sol.y(:,end);
    p1s = [p1s, tmp_sol.y(1:2:19,end)];
    p2s = [p2s, tmp_sol.y(2:2:20,end)];
    ts = [ts, ts(end)+dt];
    progressbar(i/1000)
end
figure(15201)
clf
for i = 1:10
    plot(p1s(i,:), p2s(i,:))
    hold on
end
xlabel('p1')
ylabel('p2')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Full size swarm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k_particles = 1000;
[ode_swarm_fun, q_n] = genODEswarm(sfdyn, k_particles);

%%
m = 1000;
ts = [0, zeros(1,m)];
dt = .01;
x_last = zeros(2*q_n,1);
p1s = [x_last(1:2:(q_n-1)), zeros(k_particles, m)];
p2s = [x_last(2:2:q_n),     zeros(k_particles, m)];
us = rand(q_n, m)-.5;
for i = 1:m
    tmp_sol = ode45(@(t,X) ode_swarm_fun(X(1:q_n),X((q_n+1):end),us(:,i)), [0, dt], x_last);
    x_last = tmp_sol.y(:,end);
    p1s(:,i+1) = tmp_sol.y(1:2:(q_n-1),end);
    p2s(:,i+1) = tmp_sol.y(2:2:q_n,end);
    ts(i+1) = ts(i)+dt;
    progressbar(i/m)
end
%%
figure(15202)
clf
for i = 1:k_particles
    plot(p1s(i,:), p2s(i,:))
    hold on
end
xlabel('p1')
ylabel('p2')

figure(15203)
clf
for k=fliplr([101, 201, 301, 401, 501, 601, 701, 801, 901, 1000]);
    c = (k-101)/899;
    plot(p1s(:,k), p2s(:,k),'.','Color',[c,0,1-c])
    hold on
end
xlabel('p1')
ylabel('p2')

%% Look at the bins the points fall into
[TH1, TH2] = meshgrid(linspace(-pi/2, pi/2, 20), linspace(-2*pi, 2*pi, 20));
p = zeros(size(TH1,1), size(TH1,2), m);
foo = arrayfun(@(bar) bar < TH1 & TH1 <= bar, p1s(:,end), 'UniformOutput', false)