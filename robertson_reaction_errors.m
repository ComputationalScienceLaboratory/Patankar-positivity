close all; clc

% Initial condition for Robertson's reaction
x0 = [1; 0; 0];

% Start time
t0=0;

% End time
tf=0.3;

% Array of step sizes
H = logspace(-4,-2,24);

% Error array initialization
error=zeros(length(H),1);

options = odeset('RelTol',1.e-12,'AbsTol',1.e-12, 'Jacobian', @jac_robertson);
[t,x] = ode15s(@(t,x)robertson_reaction_2(t,x),[t0,tf],x0,options);

Kmatrix_robertson = @(Y, t) calculateKmatrix(Y);

for jstep=1:length(H)
    h = H(jstep);
    
    % [t,y] = SDIRK_general(t0, tf, h, x0, @(t,x)robertson_reaction_2(t,x), 1, Kmatrix_robertson);
    [t,y] = SDIRK_general_corrected(t0, tf, h, x0, @(t,x)robertson_reaction_2(t,x), 1, Kmatrix_robertson);
    % [t,y] = RK_general(t0, tf, h, x0, @(t,y)robertson_reaction_2(t,y), 3);
    % [t,y] = RK_general_corrected(t0, tf, h, x0, @(t,x)robertson_reaction_2(t,x), 1, Kmatrix_robertson);

    error(jstep) = norm(y(:,end)' - x(end,:));
end


err = polyfit(log(H),log(error),1);

loglog(H,error);

fprintf('The empirical order is approximately: %f\n', err(1))



function dkY = jac_robertson(t, Y)  
    dkY_dY2 = [0 0 0;
        0 -3*(10^7) 0;
        0 3*(10^7) 0];

    dkY_dY3 = [0 10^4 0;
        0 -10^4 0;
        0 0 0];
    
    dkY = dkY_dY2*Y(2) + dkY_dY3*Y(3);
end

function kY = calculateKmatrix(Y, t)
    kY = [-0.04 (10^4)*Y(3) 0;
        0.04 -3*(10^7)*Y(2)-(10^4)*Y(3) 0;
        0 3*(10^7)*Y(2) 0];
end




