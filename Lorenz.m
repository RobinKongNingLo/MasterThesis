function dxdt = Lorenz(t, x)
    dxdt = zeros(3,1);
    dxdt(1) = -8/3*x(1) + x(2)*x(3);
    dxdt(2) = 1.5*(x(3)-x(2));
    dxdt(3) =  -x(1)*x(2)+28*x(2)-x(3);    
end