function example_3
% Example 3 from page 5 in "Differential Forms" by Kurt Bryan
% Example of int(w) over M = int(w(X'(t))dt) over a to b
n = 10;
t = linspace(0,2,n);

% 1-manifold, basically a line
X = [3*t;...
     t.^2;...
     5-t];

% 1-form, takes in input vector
w = @(dX) 2*X(2,:).*dX(1,:) - X(1,:).*X(3,:).*dX(2,:) + dX(3,:);


figure(13301)
clf
plot3(X(1,:), X(2,:), X(3,:))
hold on
quiver3(X(1,:), X(2,:), X(3,:), ...
        3*ones(1,n), 2*t, -ones(1,n))
quiver3(X(1,:), X(2,:), X(3,:), ...
        2*X(2,:), -X(1,:).*X(3,:), ones(1,n))
xlabel('x')
ylabel('y')
zlabel('z')

figure(13302)
clf
plot3(X(1,:), X(2,:), X(3,:))
hold on
quiver3(X(1,:), X(2,:), X(3,:), ...
        2*X(2,:), -X(1,:).*X(3,:), ones(1,n))
title('Integral over manifold')
xlabel('x')
ylabel('y')
zlabel('z')
legend('M', '1-form')
axis equal

figure(13303)
clf
plot(t, 2*X(2,:).*3 - X(1,:).*X(3,:).*2.*t -1)
title('Integral over parameter "t"')
xlabel('t')
ylabel('w')
end