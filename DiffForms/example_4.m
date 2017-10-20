function example_4
% Example 4 from page 5 in "Differential Forms" by Kurt Bryan
% Example of int(w) over M = int(w(X'(t))dt) over a to b

% Some helper functions
% Grabs the appropriate row of a vector-valued function
rSelect = @(V, IDX) V(IDX,:);

% Vector field F in R2, example
F = @(x) [x(2,:);...
          -x(1,:) - x(2,:)];
      
% 1-form, generic version
w = @(x,dx) rSelect(F(x),1).*dx(1,:) + rSelect(F(x),2).*dx(2,:);

% Manifold M in R2, parameterized by t, example
n = 10;
t = linspace(0,1,n);
M = @(x) [x.^2; x];
dM= @(x) [2*x; ones(1,length(x))];
        


%%% PLOTTING %%%
figure(14301)
clf
plot(rSelect(M(t),1), rSelect(M(t),2))
hold on
quiver(rSelect(M(t),1), rSelect(M(t),2), ...
       rSelect(F(M(t)),1), rSelect(F(M(t)),2) )
quiver(rSelect(M(t),1), rSelect(M(t),2), ...
       rSelect(dM(t),1), rSelect(dM(t),2))
xlabel('x')
ylabel('y')
zlabel('z')
legend('M', 'F', 'dM', 'Location', 'Northwest')

figure(14302)
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

figure(14303)
clf
plot(t, w(M(t),dM(t)))
title('Integral over parameter "t"')
xlabel('t')
ylabel('w')

disp('Answer is:')
disp(integral(@(x) w(M(x),dM(x)), t(1), t(end)))
end