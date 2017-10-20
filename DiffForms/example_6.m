function example_6
% Example 6 from page 6 in "Differential Forms" by Kurt Bryan
% Messing with 2-forms and 2-manifolds

% Some helper functions
% Grabs the appropriate row of a vector-valued function
rSelect = @(V, IDX) V(IDX,:);

% Vector field F in R2, example
F = @(x) [x(2,:);...
          -x(1,:) - x(3,:);...
          x(2,:)];
      
% 2-form, for curl-like stuff
w = @(x,dx) rSelect(F(x),1).*dx(1,:) + rSelect(F(x),2).*dx(2,:);

% 2-Manifold X in R3, parameterized by (u1,u2), example
n = 10;
b1 = linspace(0,1,n);
b2 = linspace(0,1,n);
[b1,b2] = meshgrid(b1,b2);
b1 = reshape(b1,[1,size(b1,1)*size(b1,2)]);
b2 = reshape(b2,[1,size(b2,1)*size(b2,2)]);
X = @(u1,u2) [u1; u1-u2; 1 - u1 + u1.*u2];
dXu1 = @(u1,u2) [ones(1, size(u1,2)); ones(1,size(u1,2)); -1+u2];
dXu2 = @(u1,u2) [zeros(1, size(u1,2)); -ones(1,size(u1,2)); u1];


%%% PLOTTING %%%
figure(15301)
clf
plot3(rSelect(X(b1,b2),1), rSelect(X(b1,b2),2), rSelect(X(b1,b2),3), '.')
hold on
quiver3(rSelect(X(b1,b2),1), rSelect(X(b1,b2),2), rSelect(X(b1,b2),3),...
       rSelect(dXu1(b1,b2),1), rSelect(dXu1(b1,b2),2), rSelect(dXu1(b1,b2),3))
quiver3(rSelect(X(b1,b2),1), rSelect(X(b1,b2),2), rSelect(X(b1,b2),3),...
       rSelect(dXu2(b1,b2),1), rSelect(dXu2(b1,b2),2), rSelect(dXu2(b1,b2),3))
quiver3(rSelect(X(b1,b2),1), rSelect(X(b1,b2),2), rSelect(X(b1,b2),3),...
       rSelect(F(X(b1,b2)),1), rSelect(F(X(b1,b2)),2), rSelect(F(X(b1,b2)),3) )
xlabel('x')
ylabel('y')
zlabel('z')
legend('M','dMu1', 'dMu2', 'F', 'Location', 'Northwest')

figure(15302)
clf
plot(rSelect(M(t),1), rSelect(M(t),2))
hold on
quiver(rSelect(M(t),1), rSelect(M(t),2), ...
       rSelect(F(M(t)),1), rSelect(F(M(t)),2) )
title('Integral over manifold')
xlabel('x')
ylabel('y')
legend('M', '1-form', 'Location', 'Northwest')
axis equal

figure(15303)
clf
plot(t, w(M(t),dM(t)))
title('Integral over parameter "t"')
xlabel('t')
ylabel('w')

disp('Answer is:')
disp(integral(@(x) w(M(x),dM(x)), t(1), t(end)))
end