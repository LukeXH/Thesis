function dynamics_integration
%DYNAMICS_INTEGRATION Just gonna take the dynamics of the smd system and
%integrate their vectors into a manifold
%   Detailed explanation goes here
m = 1;      % Mass
c = 1;      % Damping
k = 1;      % Spring Constant
w = .731/(2*pi);  % Frequency
w = 1/(2*pi);

g = 0.02;    % Input gain

fcn_phys = @(q,u) [q(2);...
                   (u - c*q(2) - k*q(1))/m];

end

