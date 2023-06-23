function dx = func(t,x)
    a = Kmatrix_func(x);
    
    dx = a*x;
end

function kY = Kmatrix_func(Y)
    kY = [-1];
end

