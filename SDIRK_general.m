function [t,y] = SDIRK_general(t0, tf, h, y0, ode_func, sdirk_method, matrix_func )
% t0 - start time
% tf - end time
% h - step size
% y0 - initial value
% ode_func - ODE function function handle
% sdirk_method - function that returns constant coefficients associated with method

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
opt.StepTolerance = 1e-10;
opt.FunctionTolerance = 1e-10;
opt.OptimalityTolerance = 1e-10;
opt.MaxIterations = 1000;
opt.MaxFunctionEvaluations = 1000;
opt.FiniteDifferenceType = 'central';

% Loop over time
for i = 2:length(t)
    dt = t(i) - t(i - 1);

    % Loop over stages
    for istage = 1:s
        U = y(:,i-1);

        for j = 1:(istage - 1)
            U = U + dt*A(istage, j)*F(:,j);
        end

        % Y(:,istage) = U + dt*gamma*ode_func(t(i), Y(:,istage));
% 		func_y = @(ys) U - ys + dt*gamma*ode_func(t(i-1) + dt*c(istage), ys);

        solver_func = @(y_stage)nonlinear_system_solver(U, dt, gamma, y_stage, matrix_func, t(i-1) + dt*c(istage));
        [y_stage, ~, exitflag, ~] = fsolve(solver_func, y(:,i-1), opt);
		if exitflag < 1
			disp 'iter';
		end
        Y(:,istage) = y_stage;

        F(:,istage) = ode_func(t(i-1) + dt*c(istage), Y(:,istage));
	end


    % Update solution
	if ~stiffaccurate
    	y_new = y(:,i-1);
    	for idx = 1:s
        	y_new = y_new + dt*b(idx)*F(:,idx);
    	end
	
    	y(:,i) = y_new;

	else
	 	y(:,i) = y_stage;

	end
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
%         gamma = .2113248654051871177454256097490213d0;
% 
%         A(1,1) = .2113248654051871177454256097490213d0;
%         A(2,1) = .2113248654051871177454256097490213d0;
%         A(2,2) = .2113248654051871177454256097490213d0;
%         A(3,1) = .2113248654051871177454256097490213d0;
%         A(3,2) = .5773502691896257645091487805019573d0;
%         A(3,3) = .2113248654051871177454256097490213d0;
% 
%         b(1)   = .2113248654051871177454256097490213d0;
%         b(2)   = .5773502691896257645091487805019573d0;
%         b(3)   = .2113248654051871177454256097490213d0;

% 		gamma = 1/4;
% 	
%         A(1,1) = 1/4;
%         A(2,1) = 1/12;
%         A(2,2) = 1/4;
%         A(3,1) = 0;
%         A(3,2) = 3/4;
%         A(3,3) = 1/4;
% 
%         b(1)   = 0;
%         b(2)   = 3/4;
%         b(3)   = 1/4;

		gamma = 1/3;
	
        A = [[1/3,0,0];[1/6,1/3,0];[5/6,-5/12,1/3]];

		b = [6/5, -1, 4/5];
    case(4)
        % 5 stages
%         gamma = .2666666666666666666666666666666667d0;
% 
%         A(1,1) = .2666666666666666666666666666666667d0;
%         A(2,1) = .5000000000000000000000000000000000d0;
%         A(2,2) = .2666666666666666666666666666666667d0;
%         A(3,1) = .3541539528432732316227461858529820d0;
%         A(3,2) = -.5415395284327323162274618585298197d-1;
%         A(3,3) = .2666666666666666666666666666666667d0;
%         A(4,1) = .8515494131138652076337791881433756d-1;
%         A(4,2) = -.6484332287891555171683963466229754d-1;
%         A(4,3) = .7915325296404206392428857585141242d-1;
%         A(4,4) = .2666666666666666666666666666666667d0;
%         A(5,1) = 2.100115700566932777970612055999074d0;
%         A(5,2) = -.7677800284445976813343102185062276d0;
%         A(5,3) = 2.399816361080026398094746205273880d0;
%         A(5,4) = -2.998818699869028161397714709433394d0;
%         A(5,5) = .2666666666666666666666666666666667d0;
% 
%         b(1)   = 2.100115700566932777970612055999074d0;
%         b(2)   = -.7677800284445976813343102185062276d0;
%         b(3)   = 2.399816361080026398094746205273880d0;
%         b(4)   = -2.998818699869028161397714709433394d0;
%         b(5)   = .2666666666666666666666666666666667d0;

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
