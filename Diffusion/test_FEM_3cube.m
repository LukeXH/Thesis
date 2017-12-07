% Script to look at how a FEM basis in 3 dimensions work. Using a 3-cube,
% assuming we are working in R3
clear

syms x y z w
c = sym('c',[8,1]); % Function coefficients
q = [1, x, y, z, x*y, x*z, y*z, x*y*z];


% Want function of the form f = sum(phi_i*w_i), where phi are basis functions and w
% are weights

% Here is the normal linear interpolation function, with unknown weights, c
f_prime = q*c;

% Now, using the (x,y,z) coordinates for each of the corners of the we get
% the weights.  The order of corners is 000,100,110,010, 001,101,111,011
B = [1, zeros(1,7);...
     1,1, zeros(1,6);...
     1,1,1,0,1, zeros(1,3);...
     1,0,1, zeros(1,5);...
     1,0,0,1, zeros(1,4);...
     1,1,0,1,0,1,0,0;...
     ones(1,8);...
     1,0,1,1,0,0,1,0];

w_prime = B*c;

% But I want the weights, w, to be free for my choosing, so lets make c a
% function of them
c_fun = inv(B)*w
     
% Now, truth be told, phi = q*D, and we want to find what D is, and since
% we know that q*c = f = phi*w, given that all are vectors in the correct
% orientation.  Given what we know, this equallity now becomes:
%           q*c = phi*w
%         q*inv(B)*w = q*D*w
% And then we have it, phi actually is phi = q*inv(B);
phi = q*inv(B);

% And the linear interpolation for a 3 cube becomes
f = phi*sym('w',[8,1])
% Note that this a Lagrange set of basis functions cuz the sum of all
% basis, phi, is equal to one: 1 = sum(phi*ones(8,1))