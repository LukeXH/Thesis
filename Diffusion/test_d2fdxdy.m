% Script looks at derivative over dxdy
clear

f = sym('f',[3,3]);
e = sym('e',[4,1]);

% e is the difference between elements:
%       e(1) = x(i+1) - x(i)
%       e(2) = x(i) - x(i-1)
%       e(3) = y(i+1) - y(i)
%       e(4) = y(i) - y(i-1)
d2fdxdy = ( e(4)       * (e(2)*f(3,3) + (e(1)-e(2))*f(3,2) - e(1)*f(3,1)) +...
           (e(3)-e(4)) * (e(2)*f(2,3) + (e(1)-e(2))*f(2,2) - e(1)*f(2,1)) -...
            e(3)       * (e(2)*f(1,3) + (e(1)-e(2))*f(1,2) - e(1)*f(1,1)) ...
           )/(4*e(1)*e(2)*e(3)*e(4));
       
% Create the convolutional matrix
c_cell = arrayfun(@(dv) diff(d2fdxdy,dv), f, 'UniformOutput', false);
c_sym = arrayfun(@(foo) foo{1}, c_cell);
pretty(c_sym)

%%% Masking
L = sym('L',[2,2]);
A = sym('A',[2,2]);
R = sym('R',[2,2]);

% pretty(L*A*R)

% Now take that convolutional matrix and apply it
p = zeros(11,11);
p(6,6) = 1;
% Assume even spacing, and flip to work with conv2
c_conv = rot90(eval(subs(c_sym, e, ones(4,1))),2);
d2pdxdy = conv2(c_conv,p);

figure(32000)
clf
p_now = p;
for i = 1:100
   surf(p_now)
   drawnow
   pause(.1)
   d2pdxdy = conv2(c_conv,p_now);
   p_now = p_now + .1*d2pdxdy(2:end-1, 2:end-1);
end