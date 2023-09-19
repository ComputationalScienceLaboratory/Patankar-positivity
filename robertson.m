%%
Ode_Function        = @(t,x)robertson_reaction(t,x);
Ode_Jacobian        = @jac_robertson;
Time_Interval       = [ 0 1.e3 ];
Y0                  = [ 1; 0; 0 ];

tol = [1.e-8 1.e-7 1.e-6 1.e-5 1.e-4];
error=zeros(length(tol),1);

integrator = matlode.rk.dirk.SDIRK_4_3_5;
options.Jacobian = Ode_Jacobian;

opts = odeset('RelTol',1.e-12,'AbsTol',1.e-12, 'Jacobian', @jac_robertson, 'NonNegative', ones(length(Y0),1));
[t,x] = ode15s(@(t,x)robertson_reaction(t,x),Time_Interval,Y0,opts);

for step=1:length(tol)
    options.AbsTol = tol(step);
    options.RelTol = tol(step);

    sol = integrator.integrate(Ode_Function, Time_Interval, Y0, options);

    error(step) = norm(sol.y(:,end)' - x(end,:));
end

err = polyfit(log(tol),log(error),1);

loglog(tol, error);
fprintf('The empirical order is approximately: %f\n', err(1))



%%
function dkY = jac_robertson(t, Y)  
    dkY_dY2 = [0 0 0;
        0 -3*(10^7) 0;
        0 3*(10^7) 0];

    dkY_dY3 = [0 10^4 0;
        0 -10^4 0;
        0 0 0];
    
    dkY = dkY_dY2*Y(2) + dkY_dY3*Y(3);
end

