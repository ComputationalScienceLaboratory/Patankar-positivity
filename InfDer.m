classdef InfDer
	%ALLENCAHN Summary of this class goes here
	%   Detailed explanation goes here
	
	properties
		full_dim
		y_dim
		y_exact
		z_exact
		u_exact
	end
	
	methods
		%NOTE: Assuming Homogenous Neauman Boundries
		function obj = InfDer()
			obj.full_dim = 2;
			obj.y_dim = 1;
			obj.y_exact = @(t) sqrt(tanh(t)+1);
			obj.z_exact = @(t) sqrt(sech(t));
			obj.u_exact = @(t) [obj.y_exact(t); obj.z_exact(t)];
			
		end
		
		function up = f_full(obj,t,u)
			up = [u(2)^4 / (2*u(1)); (-sinh(t)/(2*tanh(t) + 2)) * u(2)^3 * u(1)^2];
		end

		function J = jac(obj, t, u)
			tval = (-sinh(t)/(tanh(t) + 1));
			J = [[-u(2)^4 / (2*u(1)^2), 2*u(2)^3 / (u(1))];...
				[tval * u(2)^3*u(1),(3/2)*tval * u(2)^2*u(1)^2]];
		end

		function yp = f(obj,t,u)
			yp = u(2)^4 / (2*u(1));
		end

		function zp = g(obj,t,u)
			zp = (-sinh(t)/(2*tanh(t) + 2)) * u(2)^3 * u(1)^2;
		end
	end
end

