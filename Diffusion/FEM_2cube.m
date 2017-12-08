classdef FEM_2cube < handle
    %FEM_2CUBE Finite Element Method basis function class
    %   Ideally, this class allows for the generation of a 
    
    properties %(Access = protected)
        X % Matrix of X values
        Y % Matrix of Y values
        f % Approximated function
        x_size 
        y_size
    end
    
    methods
        function this = FEM_2cube(x_domain, y_domain)
            %FEM_2CUBE Construct an instance of this class
            %   Kinda like meshgrid
            this.x_size = length(x_domain);
            this.y_size = length(y_domain);
            [this.X, this.Y] = meshgrid(x_domain, y_domain);
            % initialize the function being approximated to zero
            this.f = zeros(this.y_size, this.x_size);
        end
        
        function res = dfdx(this)
            %METHOD1 Summary of this method goes here
            %   Gives you the first order difference along the x domain
            n = this.x_size;
            a = this.X(:,2:n-1) - this.X(:,1:n-2);
            b = this.X(:,3:n) - this.X(:,2:n-1);
            res = [(this.f(:,2)-this.f(:,1))./a(:,1),...
                   (b.*this.f(:,3:n) + (a-b).*this.f(:,2:n-1) - a.*this.f(:,1:n-2)) ./ (2.*a.*b),...
                   (this.f(:,end)-this.f(:,end-1))./b(:,end)];
        end % dfdx
        
        function res = dfdy(this)
            %METHOD1 Summary of this method goes here
            %   Gives you the first order difference along the x domain
            n = this.y_size;
            a = this.Y(2:n-1,:) - this.Y(1:n-2,:);
            b = this.Y(3:n,:) - this.Y(2:n-1,:);
            res = [(this.f(2,:)-this.f(1,:))./a(1,:);...
                   (b.*this.f(3:n,:) + (a-b).*this.f(2:n-1,:) - a.*this.f(1:n-2,:)) ./ (2.*a.*b);...
                   (this.f(end,:)-this.f(end-1,:))./b(end,:)];
        end % dfdx
        
        function res = getIntegral(this)
            % Simple central value integral
            dx = this.X(1,2:this.x_size) - this.X(1,1:this.x_size-1);
            dy = this.Y(2:this.y_size,1) - this.Y(1:this.y_size-1,1);
            
            I =   (this.f(1:this.y_size-1, 1:this.x_size-1) + ...
                   this.f(1:this.y_size-1, 2:this.x_size) + ...
                   this.f(2:this.y_size,   2:this.x_size) + ...
                   this.f(2:this.y_size,   1:this.x_size-1)...
                    )/4 .* (dy.*dy);
            res = sum(sum(I));    
        end % getIntegral
        
        function this = set_f(this, f)
            % This function the function f to a matrix value
            if any(size(this.f) ~= size(f))
               error('ERROR: input f is not the same size as internal f') 
            end
            
            this.f = f;
            
        end % set_f
        
    end
end

