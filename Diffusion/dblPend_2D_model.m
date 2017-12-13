%
clear

matlab_path = pwd;
matlab_path = matlab_path(1:strfind(pwd, 'MATLAB')+5);
addpath([matlab_path,'\Handy-Tools.git\GroupTheory'])
addpath([matlab_path,'\Handy-Tools.git\SFDyn'])

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

% sol = ode45(ode_fcn, [0, 15], [0; pi/4; 0; 0]);

%%% Simple loop to grow an RRT like tree
dt = .01;
n = 12;
ll = cell(2^n, 3); % linked list, {parent id, data, path length}
% Initial state
ll{1,1} = 0;
ll{1,2} = [0,0,0,0]'; % x,y,dx,dy
ll{1,3} = 0;
j_stop = 1;
for i = 1:n
    for j = 1:j_stop
        ll{j_stop+j,1} = j;
        q = ll{j,2};
        ll{j_stop+j,2} = q + dt*[q(3:4); accel_fcn(q(1:2), q(3:4), 50*2*(rand(2,1)-.5))];
        ll{j_stop+j,3} = 1+ll{j,3};
    end
    j_stop = j_stop + 2^(i-1);
end

%%% Plot it all out
% First retrace the paths taken
figure(15200)
clf
plot_lines = cell(2^(n-1),1);
k = 1;
for l = fliplr(2^(n-1)+1:2^n)
    parent = ll{l,1};
    plot_lines{k} = [ll{l,2}];
    while parent % when it gets to the root, 0, it while stop
       plot_lines{k} = [ll{parent,2}, plot_lines{k}];
       parent = ll{parent,1};
    end
    qs = plot_lines{k};
    s  = (1:size(qs(1,:),2))-1;
    plot3(qs(1,:), qs(2,:), s,'r')
    hold on
    k = k+1;
end
xlabel('p1')
ylabel('p2')
zlabel('epoch')
title('Random Tree of 2-Pendulum Motion')

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
xs = [];
ys = [];
ts = [];
dt = .01;
x_last = zeros(40,1);
for i=1:1000
    rng = rand(20,1)-.5;
    tmp_sol = ode45(@(t,X) multi_ode_fcn(t,X,rng), [0, dt], x_last);
    x_last = tmp_sol.y(:,end);
    xs = [xs, tmp_sol.y(1:2:19,[1,end])];
    ys = [ys, tmp_sol.y(2:2:20,[1,end])];
end
figure(15201)
clf
for i = 1:10
    plot(xs(i,:), ys(i,:))
    hold on
end