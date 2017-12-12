classdef FEM_3cube < handle
    %FEM_2CUBE Finite Element Method basis function class
    %   Ideally, this class allows for the generation of a 
    
    properties %(Access = protected)
        X % Matrix of X values
        Y % Matrix of Y values
        Z % Matrix of Z values
        f % Approximated function
        x_size 
        y_size
        z_size
    end
    
    methods
        function this = FEM_3cube(x_domain, y_domain, z_domain)
            %FEM_2CUBE Construct an instance of this class
            %   Kinda like meshgrid
            this.x_size = length(x_domain);
            this.y_size = length(y_domain);
            this.z_size = length(z_domain);
            [this.X, this.Y, this.Z] = meshgrid(x_domain, y_domain);
            % initialize the function being approximated to zero
            this.f = zeros(this.y_size, this.x_size);
        end
        
        function res = dfdx(this)
            %METHOD1 Summary of this method goes here
            %   Gives you the first order difference along the x domain
            n = this.x_size;
            a = this.X(:,3:n) - this.X(:,2:n-1);
            b = this.X(:,2:n-1) - this.X(:,1:n-2);
            res = [(this.f(:,2)-this.f(:,1))./a(:,1),...
                   (b.*this.f(:,3:n) + (a-b).*this.f(:,2:n-1) - a.*this.f(:,1:n-2)) ./ (2.*a.*b),...
                   (this.f(:,end)-this.f(:,end-1))./b(:,end)];
        end % dfdx
        
        function res = dfdy(this)
            %METHOD1 Summary of this method goes here
            %   Gives you the first order difference along the x domain
            n = this.y_size;
            a = this.Y(3:n,:) - this.Y(2:n-1,:);
            b = this.Y(2:n-1,:) - this.Y(1:n-2,:);
            res = [(this.f(2,:)-this.f(1,:))./a(1,:);...
                   (b.*this.f(3:n,:) + (a-b).*this.f(2:n-1,:) - a.*this.f(1:n-2,:)) ./ (2.*a.*b);...
                   (this.f(end,:)-this.f(end-1,:))./b(end,:)];
        end % dfdy
        
        function res = d2fdx2(this)
            %METHOD1 Summary of this method goes here
            %   Gives you the second order difference along the x domain
            %   assuming zero at boundaries
            n = this.x_size;
            a = this.X(:,3:n) - this.X(:,2:n-1);
            b = this.X(:,2:n-1) - this.X(:,1:n-2);
            res = [zeros(this.y_size,1),...
                   2*(b.*this.f(:,3:n) - (a+b).*this.f(:,2:n-1) + a.*this.f(:,1:n-2)) ./ (a.*b.*(a+b)),...
                   zeros(this.y_size,1)];
        end % d2fdx2
        
        function res = d2fdy2(this)
            %METHOD1 Summary of this method goes here
            %   Gives you the first order difference along the x domain
            n = this.y_size;
            a = this.Y(3:n,:) - this.Y(2:n-1,:);
            b = this.Y(2:n-1,:) - this.Y(1:n-2,:);
            res = [zeros(1,this.x_size);...
                   2*(b.*this.f(3:n,:) - (a+b).*this.f(2:n-1,:) + a.*this.f(1:n-2,:)) ./ (a.*b.*(a+b));...
                   zeros(1,this.x_size)];
        end % d2fdy2
        
        function res = d2fdxdy(this)
            %METHOD1 Summary of this method goes here
            %   Gives you the second order difference along the x and y domain
            %       e(1) = x(i+1) - x(i)
            %       e(2) = x(i) - x(i-1)
            %       e(3) = y(i+1) - y(i)
            %       e(4) = y(i) - y(i-1)
            e1 = [ones(this.y_size,1), this.X(:,3:end) - this.X(:,2:end-1), ones(this.y_size,1)];
            e2 = [ones(this.y_size,1), this.X(:,2:end-1) - this.X(:,1:end-2), ones(this.y_size,1)];
            e3 = [ones(1,this.x_size); this.Y(3:end,:) - this.Y(2:end-1,:); ones(1,this.x_size)];
            e4 = [ones(1,this.x_size); this.Y(2:end-1,:) - this.Y(1:end-2,:); ones(1,this.x_size)];
            f_aug = [zeros(1,this.x_size+2);...
                     zeros(this.y_size,1), this.f, zeros(this.y_size,1);...
                     zeros(1,this.x_size+2)];
            % Indicies     
            ind = {1:this.y_size, 2:this.x_size+1, 3:this.x_size+2};
            % Calculate
            res = 1./(4*e2.*e4) .*           f_aug(ind{1}, ind{1}) +...
                  (e2-e1)./(4*e1.*e2.*e4) .* f_aug(ind{1}, ind{2}) +...
                 -1./(4*e1.*e4) .*           f_aug(ind{1}, ind{3}) +...
                  (e4-e3)./(4*e2.*e3.*e4) .* f_aug(ind{2}, ind{1}) +...
                  ((e1-e2).*(e3-e4))./(4*e1.*e2.*e3.*e4) .* f_aug(ind{2}, ind{2}) +...
                  (e3-e4)./(4*e1.*e3.*e4) .* f_aug(ind{2}, ind{3}) +...
                 -1./(4*e2.*e3)           .* f_aug(ind{3}, ind{1}) +...
                  (e1-e2)./(4*e1.*e2.*e3) .* f_aug(ind{3}, ind{2}) +...
                  1./(4*e1.*e3)           .* f_aug(ind{3}, ind{3});
            
        end % d2fdxdy
        
        function res = getIntegral(this)
            % Simple central value integral
            dx = this.X(1,2:this.x_size) - this.X(1,1:this.x_size-1);
            dy = this.Y(2:this.y_size,1) - this.Y(1:this.y_size-1,1);
            
            I =   (this.f(1:this.y_size-1, 1:this.x_size-1) + ...
                   this.f(1:this.y_size-1, 2:this.x_size) + ...
                   this.f(2:this.y_size,   2:this.x_size) + ...
                   this.f(2:this.y_size,   1:this.x_size-1)...
                    )/4 .* (dy*dx);
            res = sum(sum(I));    
        end % getIntegral
        
        function this = normalize(this)
           % Takes f and normalizes it so the area under the surface is one
           % and there are no negative 
           this.f(this.f < 0) = 0;
           this.set_f(this.f/this.getIntegral());
        end
        
        function this = set_f(this, f)
            % This function the function f to a matrix value
            if any(size(this.f) ~= size(f))
               error('ERROR: input f is not the same size as internal f') 
            end
            
            this.f = f;
            
        end % set_f
        
    end
end

