function [ode_fcn, total_q_vars] = genODEswarm( sfdyn, n_particles )
%GENODESWARM Takes the Standard Form Dynamics (SFDyn class) for a single
% system and turns it into a swarm of n systems, returning a function
% handle to an ODE that can be solved with any ode solver
%   Detailed explanation goes here

accel_fcn = sfdyn.gen_accel_fcn;
m_vars    = size(sfdyn.pos_vars, 1);

q  = sym('q', [m_vars*n_particles, 1]);
dq = sym('dq', [m_vars*n_particles, 1]);
u  = sym('u', [m_vars*n_particles, 1]);

% accel_fun(q, dq, u)

f = [];
for i=1:n_particles
   ind = (m_vars*(i-1)+1):(m_vars*(i));
   f = [f; accel_fcn(q(ind), dq(ind), u(ind))];
end



ode_fcn = matlabFunction( [dq; f], 'Vars', {q,dq,u});
total_q_vars = m_vars*n_particles;
end

