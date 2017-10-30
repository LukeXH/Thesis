function varargout = ballistic2D(xy, dxy, output_type, options)
%BALLISTIC2D: general function for producing various outputs for 2D
%ballistic trajectories
%   Detailed explanation goes here

% Set parameters
if nargin == 4 && isstruct(options)
    g = options.g;
else
    g = 9.807; % m/sec^2, acceleration due to gravity
end

if nargin == 2
    output_type = 'Parabola Coefficients';
end

switch output_type
    case 'Parabola Coefficients'
        % Calculate a,b,c for y = ax^2 + bx + c
        coeff.a = -g/(2*dxy(1).^2);
        coeff.b = g * xy(1) / dxy(1).^2 + dxy(2)/dxy(1);
        coeff.c = xy(2) - g*xy(1)/(2*dxy(1).^2) - xy(1)*dxy(2)/dxy(1);
        
        varargout{1} = coeff;
        return
    
    case 'Zero Crossing'
        % Find the value of x when y = 0
        t = ( -dxy(2) - sqrt( dxy(2)^2 - 4*(-g/2)*xy(2)) )/(-g);
        varargout{1} = xy(1) + dxy(1)*t;
    
    case 'Velocity Matrix'
        % Return the 2x3 matrix that when multiplied by vector 
        % representation of position, [x,y,1]^T, gives the velocity vector
        varargout{1} = [0,        0, dxy(1);...
                       -g/dxy(1), 0, g*xy(1)/dxy(1) + dxy(2)]; 
    
    otherwise
        
end %End switch

end %End ballistic2D

