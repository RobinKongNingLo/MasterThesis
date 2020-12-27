function KDS = KDS_func(i, j, t_step, a, b, a0)
    if j <= i
        KDS = (1/(b-a)) * (1 - a0*(j*t_step - a));
    else
        KDS = (1/(b-a)) * (1 + a0*(b - j*t_step));
    end
end
