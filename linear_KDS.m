a0 = 2;
a = 0;
b = 6;
time_step = 0.01;
tspan = [a:time_step:b];
x0 = [1];
[t,x] = ode45(@(t,x) ODE4(t,x,a0), tspan, x0);
noise_snr = 20;
x = x(1,1:end-1)';
N = size(x, 2);
y = awgn(x, noise_snr, 'measured');

KV = zeros(N);
%t = (i-1)*t_step, tau = (j-1)*t_step
for i = 1:N
    for j = 1:N
        if j <= i
            KV(i,j) = (a - tspan(j))/(b-a);
        else
            KV(i,j) = (b - tspan(j))/(b-a);
        end
    end
end

C = 0;
B = 0;

for i = 1:N
    C_tmp = 0;
    for j = 1:N
        C_tmp = C_tmp + (KV(i,j)*y(j)*time_step);
    end
    C = C + C_tmp^2/N;
end
%C = C/N;

for i = 1:N
    for j = 1:N
        B = B - (1/N) * KV(i,j)*y(j)*y(i)*time_step + (1/(2*N)) * (KV(i,j)*y(j)*time_step)^2;
    end    
end

A = -B/C;


%%

x_pred = zeros(size(y));
for i = 1:N
    for j = 1:N
        x_pred(i) = x_pred(i) + KDS_func(i, j, time_step, a, b, A)*y(j)*time_step;
    end
end
figure(2)
plot(x)
hold on
plot(x_pred)