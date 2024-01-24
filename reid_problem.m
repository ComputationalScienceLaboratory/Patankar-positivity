close all; clc; clear

% Start time
t0=-2;

% End time
tf=2;

% Array of step sizes
H = logspace(-2,0,8);

error=zeros(length(H),1);

prob = InfDer();
initial_condition = prob.u_exact(t0);

options = odeset('RelTol',1.e-12,'AbsTol',1.e-12, 'Jacobian', @prob.jac, 'NonNegative', ones(length(initial_condition),1));
[t,x] = ode23t(@(t,x)prob.f_full(t,x),[t0,tf],initial_condition,options);

Kmatrix = @(Y, t) calculateKmatrix(Y,t);

for jstep=1:length(H)
    h = H(jstep);

    % [t,y] = SDIRK_general(t0, tf, h, initial_condition, @(t,x)prob.f_full(t,x), 6, @(t,y)prob.jac(t,y));
    % [t,y, ~] = SDIRK_general_clipped(t0, tf, h, initial_condition, @(t,x)prob.f_full(t,x), 3, @(t,y)prob.jac(t,y));
    % [t,y] = SDIRK_general_corrected(t0, tf, h, initial_condition, @(t,x)prob.f_full(t,x),4, Kmatrix, @(t,y)prob.jac(t,y));
    [t,y] = SDIRK_general_corrected_clipped(t0, tf, h, initial_condition, @(t,x)prob.f_full(t,x), 6, Kmatrix, @(t,y)prob.jac(t,y));

    error(jstep) = norm(y(:,end)' - x(end,:));
end

err = polyfit(log(H),log(error),1);

loglog(H,error);

fprintf('The empirical order is approximately: %f\n', err(1))


%
function kY = calculateKmatrix(Y, t)
kY = [0 Y(2)^3/(2*Y(1));
    (-sinh(t)/(2*tanh(t)+2))*Y(1)*Y(2)^3 0];
end



