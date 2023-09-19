%%
% t0 - start time
% tf - end time
% h - step size
% y0 - initial value
% ode_func - ODE function function handle
% sdirk_method - function that returns constant coefficients associated with method

%%
function [t,y] = SDIRK_general_corrected_clipped(t0, tf, h, y0, ode_func, sdirk_method, matrix_func, jacobian )

% Create time grid
Nsteps = ceil( (tf-t0) / h );
t=linspace(t0,tf,Nsteps+1);

% Initialize solution matrix
y = zeros(length(y0),Nsteps+1);
y(:,1) = y0;

% Coefficients for SDIRK method
[A, b, gamma] = SDIRK_coefficients(sdirk_method);
c = sum(A, 2);

stiffaccurate = all(b == A(end,:));

% Number of stages
s = length(b);

% Initialize stages
F = zeros(length(y0), s);
Y = zeros(length(y0), s);

opt.Display = 'off';
opt.StepTolerance = 1e-8;
opt.FunctionTolerance = 1e-8;

% Loop over time
for i = 2:length(t)
    dt = t(i) - t(i - 1);

    % Loop over stages
    for istage = 1:s
        U = y(:,i-1);

        for j = 1:(istage - 1)
            U = U + dt*A(istage, j)*F(:,j);
        end

        func_y = @(ys) U - ys + dt*gamma*ode_func(t(i-1) + dt*c(istage), ys);

        % solver_func = @(y_stage)nonlinear_system_solver(U, dt, gamma, y_stage, matrix_func, t(i-1) + dt*c(istage));

        y_stage = newton_iteration(func_y, @(Y)-eye(length(U))+ dt*gamma*jacobian(t(i-1)+dt*c(istage), Y), y(:,i-1));

        Y(:,istage) = max(y_stage,0);

        F(:,istage) = ode_func(t(i-1) + dt*c(istage), Y(:,istage));
    end

    % Update solution
    if ~stiffaccurate
        y_new = y(:,i-1) + dt*( F(:,:)*b(:) );
    else
        y_new = y_stage;
    end

    y_new = max(y_new, 0);

    for idx = 1:s
        Fmean = b(idx)*matrix_func(Y(:,idx), t(i-1) + dt*c(idx))*diag(Y(:,idx)./y_new);
    end

    y(:,i) = (eye(length(Fmean)) - dt*Fmean)\y(:,i-1);
end

end

%%
function u_f = newton_iteration(func, jac, u_0)
i = 0;
while i < 100
    c_i = -jac(u_0) \ func(u_0);
    u_f = c_i + u_0;

    i = i + 1;

    if norm(c_i) < 1.e-12
        return;
    end
    u_0 = u_f;
end
end

function Y = nonlinear_system_solver(U, dt, gamma, Y_stage, matrix_func, t)
KY = matrix_func(Y_stage, t);
Y = U + (dt*gamma*KY - eye(length(KY)))*Y_stage;
end


function [A, b, gamma] = SDIRK_coefficients(method)
switch (method)
    case(1)
        % 2 stages
        gamma = .2928932188134524755991556378951510d0;

        A(1,1) = .2928932188134524755991556378951510d0;
        A(2,1) = .7071067811865475244008443621048490d0;
        A(2,2) = .2928932188134524755991556378951510d0;

        b(1)   = .7071067811865475244008443621048490d0;
        b(2)   = .2928932188134524755991556378951510d0;
    case(2)
        % 2 stages
        x = 1 - sqrt(2)/2;
        gamma  = x;

        A = [[x,0];[1-x,x]];

        b = A(end,:);
    case(3)
        % 3 stages
        gamma = 1/3;

        A = [[1/3,0,0];[1/6,1/3,0];[5/6,-5/12,1/3]];

        b = [6/5, -1, 4/5];
    case(4)
        % 5 stages
        gamma = 1/4;

        A = [[1/4, 0, 0, 0, 0];...
            [13/20, 1/4, 0, 0, 0];...
            [580/1287, -175/5148, 1/4, 0, 0];...
            [12698/37375, -201/2990, 891/11500, 1/4, 0];...
            [944/1365, -400/819, 99/35, -575/252, 1/4]];
        b = A(end,:);

    case(5)
        % 5 stages
        gamma = .25d0;

        A(1,1) = 0.25d0;
        A(2,1) = 0.5d00;
        A(2,2) = 0.25d0;
        A(3,1) = 0.34d0;
        A(3,2) =-0.40d-1;
        A(3,3) = 0.25d0;
        A(4,1) = 0.2727941176470588235294117647058824d0;
        A(4,2) =-0.5036764705882352941176470588235294d-1;
        A(4,3) = 0.2757352941176470588235294117647059d-1;
        A(4,4) = 0.25d0;
        A(5,1) = 1.041666666666666666666666666666667d0;
        A(5,2) =-1.020833333333333333333333333333333d0;
        A(5,3) = 7.812500000000000000000000000000000d0;
        A(5,4) =-7.083333333333333333333333333333333d0;
        A(5,5) = 0.25d0;

        b(1)   =  1.041666666666666666666666666666667d0;
        b(2)   = -1.020833333333333333333333333333333d0;
        b(3)   =  7.812500000000000000000000000000000d0;
        b(4)   = -7.083333333333333333333333333333333d0;
        b(5)   =  0.250000000000000000000000000000000d0;

    case(6)
        gamma = 1/4;

        A = [[1/4, 0, 0, 0, 0];...
            [1/2, 1/4, 0, 0, 0];...
            [17/50, -1/25, 1/4, 0, 0];...
            [371/1360, -137/2720, 15/544, 1/4, 0];...
            [25/24, -49/48, 125/16, -85/12, 1/4]];
        b = A(end,:);
    otherwise
        % 2 stages
        gamma = .2928932188134524755991556378951510d0;

        A(1,1) = .2928932188134524755991556378951510d0;
        A(2,1) = .7071067811865475244008443621048490d0;
        A(2,2) = .2928932188134524755991556378951510d0;

        b(1)   = .7071067811865475244008443621048490d0;
        b(2)   = .2928932188134524755991556378951510d0;
end

end