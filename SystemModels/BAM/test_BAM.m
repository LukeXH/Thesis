% This is an attempt to recreate the dynamics of the Tethered ball from the
% end of the bam.

clear

matlab_path = pwd;
matlab_path = matlab_path(1:strfind(pwd, 'MATLAB')+5);
addpath([matlab_path,'\Handy-Tools.git\GroupTheory'])
addpath([matlab_path,'\Handy-Tools.git\SFDyn'])
addpath([matlab_path,'\Thesis.git\SystemModels\Trajectories'])

% % q = 180 is in direction of gravity, q increments clocwise.
% % Global coordinates of ball, from pivot of throwing arm
% x = sym('x', [2,1]);
% y = sym('y', [2,1]);
% 
% %Arm coordinates, angle and angular velocity
% q = sym('q',[2,1]);
% syms I g m r b Cd

% Equations of motion for the arm.
dq = @(q,r,g,I,b) [q(2);...
                   r*g/I*sin(q(1)) - b/I*q(2)];
% Transform from arm state space coordiantes to global positional
% coordiantes
xy_arm = @(q,r) [r*sin(q(1));...
                 r*cos(q(1))];
      
% Ballistic equations for the ball, ODE
dxy_ball = @(x,y,g) [x(2);...
                     y(2);...
                     0;...
                    -g];
       
%%%
transitionEvent = @(t,X)