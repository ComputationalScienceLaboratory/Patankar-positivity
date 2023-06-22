function [t,y] = RK_general_corrected(t0, tf, h, y0, ode_func, rk_method, matrix_func)
% t0 - start time
% tf - end time
% h - step size
% y0 - initial value
% ode_func - ODE function function handle
% rk_method - function that returns constant coefficients associated with method

% Create time grid
Nsteps = ceil( (tf-t0) / h );
t=linspace(t0,tf,Nsteps+1);

% Initialize solution matrix
y = zeros(length(y0),Nsteps+1);
y(:,1) = y0;

% Coefficients for Runge-Kutta method
[A, B, C] = RK_coefficients(rk_method);

% Number of stages
s = length(B);

% Loop over time
for i = 2:length(t)
    dt = t(i) - t(i - 1);
    Y = zeros(length(y0), s);
    Y(:, 1) = y(:, i-1);

    % Loop over stages
    for istage = 2:s
        arg = t(i-1) + C(istage)*dt;
        sum = zeros(length(y0), 1);

        for j = 1:(istage - 1)
            sum = sum + A(istage, j)*ode_func(arg, Y(:,j));
        end

        Y(:,istage) = Y(:,istage - 1) + dt*sum;
    end

    % Update solution
    y_new = y(:,i-1);
    for idx = 1:s
        time = t(i-1) + C(idx)*dt;
        y_new = y_new + dt*B(idx)*ode_func(time, Y(:,idx));
    end

    

    % Correction
    for idx = 1:s
        Kmean = B(idx)*matrix_func(Y(:,idx), t(i))*diag(Y(:,idx)./y_new);
    end

    y(:,i) = (eye(length(y0)) - dt*Kmean)\y(:,i-1);
end

end


function [A, B, C] = RK_coefficients(method)
switch (method)
    case(1)
        % 2nd order
        C = [0 1/2];
        A = [0 0; 1/2 0];
        B = [0 1];

    case(2)
        % 3rd order Heun
        C = [0 1/3 2/3];
        A = [0 0 0; 1/3 0 0; 0 2/3 0];
        B = [1/4 0 3/4];

    case(3)
        % 4th order
        C = [0 1/2 1/2 1];
        A = [0 0 0 0; 1/2 0 0 0; 0 1/2 0 0; 0 0 1 0];
        B = [1/6 1/3 1/3 1/6];

    otherwise
        % 2nd order
        C = [0 1/2];
        A = [0 0; 1/2 0];
        B = [0 1];
end

end
