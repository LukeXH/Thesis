% Script to look at how a FEM basis in 2 dimensions work. Using a 3-cube,
% assuming we are working in R3
clear

syms x y w
c = sym('c',[4,1]); % Function coefficients
q = [1, x, y, x*y];


% Want function of the form f = sum(phi_i*w_i), where phi are basis functions and w
% are weights

% Here is the normal linear interpolation function, with unknown weights, c
f_prime = q*c;

% Now, using the (x,y,z) coordinates for each of the corners of the we get
% the weights.  The order of corners is 00,10,11,01
B = [1,0,0,0;...
     1,1,0,0;...
     1,1,1,1;...
     1,0,1,0];

w_prime = B*c;

% But I want the weights, w, to be free for my choosing, so lets make c a
% function of them
c_fun = inv(B)*w;
     
% Now, truth be told, phi = q*D, and we want to find what D is, and since
% we know that q*c = f = phi*w, given that all are vectors in the correct
% orientation.  Given what we know, this equallity now becomes:
%           q*c = phi*w
%         q*inv(B)*w = q*D*w
% And then we have it, phi actually is phi = q*inv(B);
phi = q*inv(B);

% And the linear interpolation for a 3 cube becomes
f = phi*sym('w',[4,1])
% Note that this a Lagrange set of basis functions cuz the sum of all
% basis, phi, is equal to one: 1 = sum(phi*ones(8,1))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now let's look at the idea of keeping this smooth, at least to first
% order, of the function between between 2-cube elements
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%First, start with a simple 3x3 points of evaluation
n = 3;
[U,V] = meshgrid(0:n-1,0:n-1);
% And then lets say that here's what the function looks like Z at the
% following points
Z = [1, 1, 1;...
     1, 2, 1;...
     1, 1, 1];

% Now let's look at just one corner of this
m = 10;
[U1,V1] = meshgrid(linspace(0,1,m), linspace(0,1,m));
W = [Z(1,1), Z(1,2), Z(2,2), Z(2,1)]';
phi_fun = matlabFunction(phi);
Z1_interp = arrayfun(@(a,b)phi_fun(a,b)*W,U1,V1)


%%% Plotting!
% This one is how straight up linear interpolation looks like
figure(23300)
clf
surf(U,V,Z)
title('Flat Linear Interp')

figure(23301)
clf
surf(U1,V1,Z1_interp)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Testing the FEM_2cube class.  ...Finally.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = FEM_2cube(linspace(0,1,21), linspace(1,2,21));
p.f(11,11) = 1;
p.set_f(p.f/p.getIntegral())

figure(23301)
clf
surf(p.X,p.Y,p.f)
for k = 1:100
%     dp = -(p.dfdx + p.dfdy) + p.d2fdx2 + p.d2fdy2 + p.d2fdxdy;
%     dp = -p.dfdy + .01*(-p.d2fdxdy + p.d2fdx2 + p.d2fdy2);
    dp = -p.dfdy + .01*(p.d2fdx2 + p.d2fdy2);
    p.set_f(p.f + dp*0.01);
    p.f(p.f < 0) = 0;
    p.set_f(p.f/p.getIntegral());
    surf(p.X,p.Y,p.f)
    drawnow
    pause(0.1)
end
xlabel('x')
ylabel('y')

