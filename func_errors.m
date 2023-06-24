close all; clc; clear

% Initial condition
x0 = [1];

% Start time
t0=0;

% End time
tf=0.2;

% Array of step sizes
H = logspace(-3,-1, 5);

% Error array initialization
error=zeros(length(H),1);

options = odeset('RelTol',1.e-12,'AbsTol',1.e-12, 'Jacobian', @jac_mapk, 'NonNegative', ones(length(x0),1));
[t,x] = ode15s(@(t,x)func(t,x),[t0,tf],x0,options);

Kmatrix_func = @(Y, t) calculateKmatrix(Y);

for jstep=1:length(H)
    h = H(jstep);
    
    [t,y] = SDIRK_general(t0, tf, h, x0, @(t,x)func(t,x), 3, Kmatrix_func);
    % [t,y] = SDIRK_general_corrected(t0, tf, h, x0, @(t,x)func(t,x), 1, Kmatrix_func);
    % [t,y] = RK_general(t0, tf, h, x0, @(t,x)func(t,x), 4);
    % [t,y] = RK_general_corrected(t0, tf, h, x0, @(t,x)func(t,x), 2, Kmatrix_func);

    error(jstep) = norm(y(:,end)' - x(end,:));
%     error(jstep) = norm(y(:,end)' - exp(-tf));
end

err = polyfit(log(H),log(error),1);

loglog(H,error);

fprintf('The empirical order is approximately: %f\n', err(1))

function dkY = jac_mapk(t, Y)    
    dkY = [0];
end


function kY = calculateKmatrix(Y, t)
    kY = [-1];
end



