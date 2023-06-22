close all; clc

% x0 = [9.906E+01, 6.624E+08, 5.326E+11, 8.725E+08, 2.240E+08];
% fix = [8.120E+16, 1.697E+16];

% Initial condition for Stratospheric reaction
x0 = [9.906E+01; 6.624E+08; 5.326E+11; 1.697E+16; 8.725E+08; 2.240E+08];


% Start time
t0=0;

% End time
tf=0.3;

% Array of step sizes
H = logspace(-3,-1,6);

% Error array initialization
error=zeros(length(H),1);

options = odeset('RelTol',1.0e-6,'AbsTol',1.0e-3, 'Jacobian', @jac_stratospheric);
[t,x] = ode15s(@(t,x)stratospheric_reaction_2(t,x),[t0, tf],x0,options);

Kmatrix_stratosperic = @(Y, t) calculateKmatrix(Y, t);

for jstep=1:length(H)
    h = H(jstep);

    [t,y] = SDIRK_general(t0, tf, h, x0, @(t,x)stratospheric_reaction_2(t,x), 3, Kmatrix_stratosperic);
    % [t,y] = SDIRK_general_corrected(t0, tf, h, x0, @(t,x)stratospheric_reaction_2(t,x), 3, Kmatrix_stratosperic);
    % [t,y] = RK_general(t0, tf, h, x0, @(t,x)stratospheric_reaction_2(t,x), 3);
    % [t,y] = RK_general_corrected(t0, tf, h, x0, @(t,x)stratospheric_reaction_2(t,x), 3, Kmatrix_stratosperic);
    
    error(jstep) = norm(y(:,end)' - x(end,:));
end

err = polyfit(log(H),log(error),1);

loglog(H,error);

fprintf('The empirical order is approximately: %f\n', err(1))


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

function kY = calculateKmatrix(Y, t)
    % K matrix for stratospheric example
    [k] = stratospheric_parameters(t);

    kY = [(-k(6)-k(7)*Y(3)) 0 k(5) 0 0 0;
        k(6) (-k(2)*Y(4)-k(4)*Y(3)-k(9)*Y(6)) k(3) 2*k(1) 0 k(10);
        0 (1/3)*k(2)*Y(4) (-k(3)-k(5)-k(4)*Y(2)-k(7)*Y(1)-k(8)*Y(5)) (2/3)*k(2)*Y(2) 0 0;
        (1/2)*k(7)*Y(3) (k(4)*Y(3)+(1/2)*k(9)*Y(6)) k(3)+k(5)+k(4)*Y(2)+k(7)*Y(1)+k(8)*Y(5)+(1/2)*k(7)*Y(1) (-k(1)-k(2)*Y(2)) 0 (1/2)*k(9)*Y(2);
        0 0 0 0 (-k(8)*Y(3)) k(10)+k(9)*Y(2);
        0 0 0 0 k(8)*Y(3) (-k(10)-k(9)*Y(2))];
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





