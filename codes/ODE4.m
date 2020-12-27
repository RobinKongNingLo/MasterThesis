function dxdt = ODE4(t, x, a0, a1, a2, a3)
    dxdt = zeros(4,1);
    dxdt(1) = x(2);
    dxdt(2) = x(3);
    dxdt(3) = x(4);
    dxdt(4) = - a3.*x(4) - a2.*x(3) - a1.*x(2) - a0.*x(1);    
end
