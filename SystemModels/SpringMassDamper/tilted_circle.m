clear

n = 100;
th = linspace(0, 2*pi, n);

a = [0,0,0]';
b = [1,0,0]';
c = [0,1,1]';
R = 0.2;

r = repmat(a, 1, 100) + b*R*cos(th) + c*R*sin(th);

figure(12399)
clf
plot3(r(1,:), r(2,:), r(3,:))
xlabel('x')
ylabel('y')
zlabel('z')