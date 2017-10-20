function dynamics_fractal
%DYNAMICS_FRACTAL Trying to recreate manifold for by looking at the horizon
%for a variety of actions.
m = 1;      % Mass
c = 1;      % Damping
k = 1;      % Spring Constant
w = .731/(2*pi);  % Frequency
w = 1/(2*pi);

g = 0.02;    % Input gain
dt = 0.01;

us = g*[-1,0,1];

fcn_phys = @(q,u) [q(2);...
                   (u - c*q(2) - k*q(1))/m];

               
%%% Simple branching fractal
figure(15301)
clf
plot3([0],[0],[0],'kx')
xlabel('u')
ylabel('x')
zlabel('dx')
hold on
%%% Plot the fractal
ends = [0;0;0];
for i = 1:6; % Fractal depth
    disp(sprintf('Iteration %f', i))
    new_ends = zeros(3, size(ends,2)*3);
    for n = 1:size(ends,2)
        for m=1:3
           ind = (3*(n-1)+m);
           new_ends(:,ind) = [us(m);...
                              dt*fcn_phys(ends(2:3,n),us(m)) + ends(2:3,n)];
           line([ends(1,n), new_ends(1,ind)],...
                [ends(2,n), new_ends(2,ind)],...
                [ends(3,n), new_ends(3,ind)]);
        end
    end
    ends = new_ends;
end

%%% "Fan" fractal
figure(15302)
clf
plot3([0],[0],[0],'kx')
xlabel('u')
ylabel('x')
zlabel('dx')
hold on
%%% Plot the fractal
wid = 6;
depth = 3;
us_mod = linspace(us(1), us(end), wid);
ends = [0;0;0];
colors = [linspace(1,0,depth)', linspace(.3,0,depth)', linspace(0,1,depth)'];
for i = 1:depth; % Fractal depth
    disp(sprintf('Iteration %f', i))
    new_ends = zeros(3, size(ends,2)*3);
    for n = 1:size(ends,2)
        for m=1:wid
           ind = (wid*(n-1)+m);
           new_ends(:,ind) = [us_mod(m);...
                              dt*fcn_phys(ends(2:3,n),us_mod(m)) + ends(2:3,n)];
          
        end
        plot3(new_ends(1,:), new_ends(2,:), new_ends(3,:),'.', 'Color', colors(i,:));
    end
    ends = new_ends;
end


end

