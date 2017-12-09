% Script looks at derivative over dxdy
clear

f = sym('f',[3,3]);
e = sym('e',[4,1]);

d2fdxdy = ( e(4)       * (e(2)*f(3,3) + (e(1)-e(2))*f(3,2) - e(1)*f(3,1)) +...
           (e(3)-e(4)) * (e(2)*f(2,3) + (e(1)-e(2))*f(2,2) - e(1)*f(2,1)) -...
            e(3)       * (e(2)*f(1,3) + (e(1)-e(2))*f(1,2) - e(1)*f(1,1)) ...
           )/(4*e(1)*e(2)*e(3)*e(4));
       
c_cell = arrayfun(@(dv) diff(d2fdxdy,dv), f, 'UniformOutput', false)


%%% Masking
L = sym('L',[2,2]);
A = sym('A',[2,2]);
R = sym('R',[2,2]);

pretty(L*A*R)