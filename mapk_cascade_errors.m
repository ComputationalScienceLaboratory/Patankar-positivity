close all; clc; clear

% Initial condition for MAPK cascade
x0 = [0.1; 0.175; 0.15; 1.15; 0.81; 0.5];

% Start time
t0=0;

% End time
tf=0.2;

% Array of step sizes
H = logspace(-3,-2,12);

% Error array initialization
error=zeros(length(H),1);

options = odeset('RelTol',1.e-12,'AbsTol',1.e-12, 'Jacobian', @jac_mapk, 'NonNegative', ones(length(x0),1));
[t,x] = ode15s(@(t,x)mapk_cascade_2(t,x),[t0,tf],x0,options);

Kmatrix_mapk = @(Y, t) calculateKmatrix(Y);

for jstep=1:length(H)
    h = H(jstep);
    
    % [t,y] = SDIRK_general(t0, tf, h, x0, @(t,x)mapk_cascade_2(t,x), 1, Kmatrix_mapk);
    % [t,y] = SDIRK_general_corrected(t0, tf, h, x0, @(t,x)mapk_cascade_2(t,x), 5, Kmatrix_mapk);
    % [t,y] = SDIRK_general_clipped(t0, tf, h, x0, @(t,x)mapk_cascade_2(t,x), 3, Kmatrix_mapk);
    % [t,y] = RK_general(t0, tf, h, x0, @(t,x)mapk_cascade_2(t,x), 3);
    % [t,y] = RK_general_clipped(t0, tf, h, x0, @(t,x)mapk_cascade_2(t,x), 2, Kmatrix_mapk);
    % [t,y] = RK_general_corrected(t0, tf, h, x0, @(t,x)mapk_cascade_2(t,x), 3, Kmatrix_mapk);

    error(jstep) = norm(y(:,end)' - x(end,:));
end

err = polyfit(log(H),log(error),1);

loglog(H,error);

fprintf('The empirical order is approximately: %f\n', err(1))

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


function kY = calculateKmatrix(Y, t)
    [alpha, k] = mapk_parameters();
    kY = [(-k(7)-k(1)*Y(2)) 0 0 k(2) 0 k(6); 
          0 (-k(1)*Y(1)) k(5) 0 0 0; 
          0 0 (-k(3)*Y(1)-k(5)) k(2) k(4) 0; 
          (1-alpha)*(k(1)*Y(2)) (alpha*k(1)*Y(1)) 0 (-k(2)) 0 0; 
          0 0 (k(3)*Y(1)) 0 (-k(4)) 0; 
          k(7) 0 0 0 0 (-k(6))];
end



