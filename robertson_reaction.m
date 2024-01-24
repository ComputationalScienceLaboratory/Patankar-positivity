function dx = robertson_reaction(t,x)
    a = Kmatrix_robertson(x);
    
    dx = a*x;
end

function kY = Kmatrix_robertson(Y)
    % K matrix for robertson example
    kY = [-0.04 (10^4)*Y(3) 0;
        0.04 -3*(10^7)*Y(2)-(10^4)*Y(3) 0;
        0 3*(10^7)*Y(2) 0];
end

