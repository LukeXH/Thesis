clear

% This script is my fiddling around with parallel transport, trying to get
% a grasp on it.  Currently I'm looking at the geodesics on a sphere.  The
% sphere's center is (1,1,1) and has a radius of 1.  We'll look at
% transport at the equatorial geodesic and then on a curve at a high
% latitude

% Defining the sphere
r = 1;          % Radius
o = [1; 1; 1];  % Origin
t = linspace(0, 2*pi, 40); % Theta parameterization
xz_circ = [1 + r*cos(t); 1*ones(1,size(t,2)); 1 + r*sin(t)];
yz_circ = [1*ones(1,size(t,2)); 1 + r*cos(t); 1 + r*sin(t)];

% Equatorial geodesic
geo = [1 + r*cos(t); 1 + r*sin(t); 1*ones(1,size(t,2))];

% curve path
psi = pi/4;
c = [o(1) + r*cos(psi)*cos(t); o(2) + r*cos(psi)*sin(t); (o(3)+r*sin(psi))*ones(1,size(t,2))];
dcdt = [-r*cos(psi)*sin(t); r*cos(psi)*cos(t); zeros(1,size(t,2))];
ddcdtt = [-r*cos(psi)*cos(t); -r*cos(psi)*sin(t); zeros(1,size(t,2))];

%Piecewise paths
th = linspace(0, pi/2, 20);
phi= linspace(0, pi/4, 20);
c2 = [o(1) + r*cos(phi)*cos(0), o(1) + r*cos(phi(end))*cos(th), o(1) + r*cos(fliplr(phi))*cos(th(end)), o(1) + r*cos(phi(1))*cos(fliplr(th));...
      o(2) + r*cos(phi)*sin(0), o(2) + r*cos(phi(end))*sin(th), o(2) + r*cos(fliplr(phi))*sin(th(end)), o(2) + r*cos(phi(1))*sin(fliplr(th));...
      o(3)+r*sin(phi),          o(3)+r*sin(phi(end))*ones(1,size(th,2)), o(3)+r*sin(fliplr(phi)),       o(3)+r*sin(phi(1))*ones(1,size(th,2))];

% Tangent Space (Asssume sphere at origin)
dp1 = @(x,y,z) [-y; x; 0.*y.*x];
dp2 = @(x,y,z) [-z.*x./sqrt(x.^2+y.^2); -z.*y./sqrt(x.^2+y.^2); sqrt(r.^2 - z.^2)];
TpM = @(x,y,z) [dp1(x,y,z); dp2(x,y,z)]; % 1-3: dp1, 4-6: dp2

% Vector field, V
V = r/4 * eye(3);

%%% Display
figure(123400)
clf
plot3(xz_circ(1,:), xz_circ(2,:), xz_circ(3,:), 'b')
hold on
plot3(yz_circ(1,:), yz_circ(2,:), yz_circ(3,:), 'b')

% The geodesic
plot3(geo(1,:), geo(2,:), geo(3,:), 'r')
T = TpM(geo(1,:)-o(1), geo(2,:)-o(2), geo(3,:)-o(3));
T = [T(1,:)./sqrt(sum(T(1:3,:).^2));...
     T(2,:)./sqrt(sum(T(1:3,:).^2));...
     T(3,:)./sqrt(sum(T(1:3,:).^2));...
     T(4,:)./sqrt(sum(T(4:6,:).^2));...
     T(5,:)./sqrt(sum(T(4:6,:).^2));...
     T(6,:)./sqrt(sum(T(4:6,:).^2))];
T1u = T(1:3,1)/norm(T(1:3,1));
T1v = T(4:6,1)/norm(T(4:6,1));
patch_verticies = [r/4*(-T1u-T1v), r/4*(T1u-T1v), r/4*(T1u+T1v), r/4*(-T1u+T1v)] + repmat(geo(:,1), [1,4]);
patch(patch_verticies(1,:), patch_verticies(2,:), patch_verticies(3,:), [0, .5, .9], 'FaceAlpha', .5)
% quiver3(geo(:,1), geo(:,2), geo(:,3), T(:,1), T(:,2), T(:,3))
% quiver3(geo(:,1), geo(:,2), geo(:,3), T(:,4), T(:,5), T(:,6))

% The curve
plot3(c(1,:), c(2,:), c(3,:), 'm')
T = TpM(c(1,:)-o(1), c(2,:)-o(2), c(3,:)-o(3));
T = [T(1,:)./sqrt(sum(T(1:3,:).^2));...
     T(2,:)./sqrt(sum(T(1:3,:).^2));...
     T(3,:)./sqrt(sum(T(1:3,:).^2));...
     T(4,:)./sqrt(sum(T(4:6,:).^2));...
     T(5,:)./sqrt(sum(T(4:6,:).^2));...
     T(6,:)./sqrt(sum(T(4:6,:).^2))];
T1u = T(1:3,1)/norm(T(1:3,1));
T1v = T(4:6,1)/norm(T(4:6,1));
patch_verticies = [r/4*(-T1u-T1v), r/4*(T1u-T1v), r/4*(T1u+T1v), r/4*(-T1u+T1v)] + repmat(c(:,1), [1,4]);
patch(patch_verticies(1,:), patch_verticies(2,:), patch_verticies(3,:), [0, .5, .9], 'FaceAlpha', .5)
% quiver3(c(1,:), c(2,:), c(3,:), T(1,:), T(2,:), T(3,:))
% quiver3(c(1,:), c(2,:), c(3,:), T(4,:), T(5,:), T(6,:))

% The curve's parameter derivitive
quiver3(c(1,1), c(2,1), c(3,1),  dcdt(1,1), dcdt(2,1), dcdt(3,1), 'k')
quiver3(c(1,1), c(2,1), c(3,1),  ddcdtt(1,1), ddcdtt(2,1), ddcdtt(3,1), 'k')

% Broken Curve/Piecewise path
% plot3(c2(1,:), c2(2,:), c2(3,:), 'Color', [0, .5, 0], 'Marker', 'o')

axis equal
axis square
title('Is this Parallel Transport?')