%%
Ode_Function        = @(t,x)stratospheric_reaction(t,x);
Ode_Jacobian        = @jac_stratospheric;
Time_Interval       = [ 12*3600 36*3600 ];
Y0                  = [ 9.906E+01; 6.624E+08; 5.326E+11; 1.697E+16; 8.725E+08; 2.240E+08 ];

tol = [1.e-8 1.e-7 1.e-6 1.e-5 1.e-4];
error=zeros(length(tol),1);

integrator = matlode.rk.dirk.SDIRK_4_3_5;
options.Jacobian = Ode_Jacobian;

opts = odeset('RelTol',1.e-12,'AbsTol',1.e-12, 'Jacobian', @jac_stratospheric, 'NonNegative', ones(length(Y0),1));
[t,x] = ode15s(@(t,x)stratospheric_reaction(t,x),Time_Interval,Y0,opts);

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
function [k] = stratospheric_parameters(t)
    Tr=4.5;
    Ts=19.5;
    Tl=mod((t/3600),24);

    if((Tr<=Tl) && (Tl<=Ts))
        Ttmp = (2*Tl-Tr-Ts)/(Ts-Tr);
        sigma=(1/2)+(1/2)*cos(pi*(abs(Ttmp))*(Ttmp));
    else
        sigma=0;
    end

    k=[(2.643E-10)*(sigma^3), 8.018E-17, (6.120E-04)*sigma, 1.576E-15, (1.070E-03)*(sigma^2), 7.110E-11, 1.200E-10, 6.062E-15, 1.069E-11, (1.289E-02)*sigma];

end

function dkY = jac_stratospheric(t, Y)
    [k] = stratospheric_parameters(t);
    dkY_dY1 = [0 0 0 0 0 0;
        0 0 0 0 0 0;
        0 0 -k(7) 0 0 0;
        0 0 (3/2)*k(7) 0 0 0;
        0 0 0 0 0 0;
        0 0 0 0 0 0];
    
    dkY_dY2 = [0 0 0 0 0 0;
        0 0 0 0 0 0;
        0 0 -k(4) (2/3)*k(2) 0 0;
        0 0 k(4) -k(2) 0 (1/2)*k(9);
        0 0 0 0 0 k(9);
        0 0 0 0 0 -k(9)];

    dkY_dY3 = [-k(7) 0 0 0 0 0;
        0 -k(4) 0 0 0 0;
        0 0 0 0 0 0;
        (1/2)*k(7) k(4) 0 0 0 0;
        0 0 0 0 -k(8) 0;
        0 0 0 0 k(8) 0];

    dkY_dY4 = [0 0 0 0 0 0;
        0 -k(2) 0 0 0 0;
        0 (1/3)*k(2) 0 0 0 0;
        0 0 0 0 0 0;
        0 0 0 0 0 0;
        0 0 0 0 0 0];

    dkY_dY5 = [0 0 0 0 0 0;
        0 0 0 0 0 0;
        0 0 -k(8) 0 0 0;
        0 0 k(8) 0 0 0;
        0 0 0 0 0 0;
        0 0 0 0 0 0];

    dkY_dY6 = [0 0 0 0 0 0;
        0 -k(9) 0 0 0 0;
        0 0 0 0 0 0;
        0 (1/2)*k(9) 0 0 0 0;
        0 0 0 0 0 0;
        0 0 0 0 0 0];
    
    dkY = dkY_dY1*Y(1) + dkY_dY2*Y(2) + dkY_dY3*Y(3) + dkY_dY4*Y(4) + dkY_dY5*Y(5) + dkY_dY6*Y(6);
end

