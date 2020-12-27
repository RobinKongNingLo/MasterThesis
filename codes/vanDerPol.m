function dxdt = vanDerPol(t, x)
    dxdt = zeros(2,1);
    %dxdt(1) = 2 * (x(1) - 1/3 * x(1)^3 - x(2));
    %dxdt(2) = 0.5 * x(1);
    dxdt(1) = x(2);
    dxdt(2) = 0.5*(1-x(1)^2)*x(2)-x(1);
end
