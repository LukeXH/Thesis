function func_handle = transitionEvents( selection, param )
%TRANSITIONEVENTS Summary of this function goes here
%   Detailed explanation goes here



function [value, isterminal, direction] = release(T, Y)
    value      = (Y(3) == 0.1);
    isterminal = 1;   % Stop the integration
    direction  = 0;
end


end

