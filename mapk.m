close all; clc; clear
%%
Ode_Function        = @(t,x)mapk_cascade(t,x);
Ode_Jacobian        = @jac_mapk;
Time_Interval       = [ 0 200 ];
Y0                  = [ 0.1; 0.175; 0.15; 1.15; 0.81; 0.5 ];

tol = [1.e-10 1.e-9 1.e-8 1.e-7 1.e-6];
error=zeros(length(tol),1);

steps=zeros(length(tol),1);

integrator = matlode.rk.dirk.SDIRK_2_1_2;
options.Jacobian = Ode_Jacobian;


sol = integrator.integrate(Ode_Function, Time_Interval, Y0, options);

for step=1:length(tol)
    % options.ErrNorm = matlode.errnorm.InfNorm(tol(step), tol(step));

    opts = odeset('RelTol',tol(step),'AbsTol',tol(step), 'Jacobian', @jac_mapk, 'NonNegative', ones(length(Y0),1));
    [t,x] = ode23t(@(t,x)mapk_cascade(t,x), Time_Interval, Y0, opts);
    

    error(step) = norm(sol.y(:,end)' - x(end,:));
    % steps(step) = sol.stats.nSteps;
end

err = polyfit(log(tol),log(error),1);

loglog(tol, error);
fprintf('The empirical order is approximately: %f\n', err(1))


%%
function [alpha,k] = mapk_parameters
alpha=1;
k=[100/3; 1/3; 50; 0.5; 10/3; 0.1; 7/10];
end

function dkY = jac_mapk(t, Y)
[alpha,k] = mapk_parameters;
dkY_dY1 = [0 0 0 0 0 0;
    0 (-k(1)) 0 0 0 0;
    0 0 (-k(3)) 0 0 0;
    0 (alpha*k(1)) 0 0 0 0;
    0 0 (k(3)) 0 0 0;
    0 0 0 0 0 0];

dkY_dY2 = [(-k(1)) 0 0 0 0 0;
    0 0 0 0 0 0;
    0 0 0 0 0 0;
    (1-alpha)*(k(1)) 0 0 0 0 0;
    0 0 0 0 0 0;
    0 0 0 0 0 0];

dkY = dkY_dY1*Y(1) + dkY_dY2*Y(2);
end

