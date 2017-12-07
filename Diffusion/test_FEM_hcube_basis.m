% This script displays a Lagrange interpolation basis for a 2D hypercube
% (a square)
clear

% Weights for each basis, the basis are centered a (0,0), (1,0), (1,1),
% (0,1)
w = [1,1,1,1];

% Basis functions
p1 = @(x,y) 1 - x - y + x.*y;
p2 = @(x,y) x - x.*y;
p3 = @(x,y) x.*y;
p4 = @(x,y) y - x.*y;

%%% Visualization
n = 5;
[X,Y] = meshgrid(linspace(0,1,n),linspace(0,1,n));

figure(111100)
clf
subplot(2,2,1)
surf(X,Y, w(1)*p1(X,Y) + w(2)*p2(X,Y) + w(3)*p3(X,Y) + w(4)*p4(X,Y))
subplot(2,2,2)
w = [2,1,1,1];
surf(X,Y, w(1)*p1(X,Y) + w(2)*p2(X,Y) + w(3)*p3(X,Y) + w(4)*p4(X,Y))
subplot(2,2,3)
w = [2,2,1,1];
surf(X,Y, w(1)*p1(X,Y) + w(2)*p2(X,Y) + w(3)*p3(X,Y) + w(4)*p4(X,Y))
subplot(2,2,4)
w = [2,2,2,1];
surf(X,Y, w(1)*p1(X,Y) + w(2)*p2(X,Y) + w(3)*p3(X,Y) + w(4)*p4(X,Y))