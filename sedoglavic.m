function dxdt = sedoglavic(t, x)
    dxdt = zeros(3,1);
    dxdt(1) = x(2) / x(1);
    dxdt(2) = x(3) / x(2);
    dxdt(3) = x(1) * 0.5;
end