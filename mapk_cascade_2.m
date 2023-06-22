function dx = mapk_cascade_2(t,x)
    a = Kmatrix_mapk(x);
    
    dx = a*x;
end

function [alpha,k] = mapk_parameters
    alpha=1;
    k=[100/3; 1/3; 50; 0.5; 10/3; 0.1; 7/10];
end

function kY = Kmatrix_mapk(Y)
    % K matrix for MAPK example
    [alpha,k] = mapk_parameters;
    kY = [(-k(7)-k(1)*Y(2)) 0 0 k(2) 0 k(6); 
        0 (-k(1)*Y(1)) k(5) 0 0 0; 
        0 0 (-k(3)*Y(1)-k(5)) k(2) k(4) 0; 
        (1-alpha)*(k(1)*Y(2)) (alpha*k(1)*Y(1)) 0 (-k(2)) 0 0; 
        0 0 (k(3)*Y(1)) 0 (-k(4)) 0; 
        k(7) 0 0 0 0 (-k(6))];
end

