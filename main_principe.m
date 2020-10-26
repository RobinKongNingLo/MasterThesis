%System initialization
m = 200;
c1 = 0.4;
c2 = 0.84;
c3 = 6.0;
c4 = 1.0;
x = zeros(2, m+2); %Data: [x0, x1, ..., xm, xm+1]
x(:,1) = [1;0]; %x0 = [1 0]^T
noise_snr = -3;
%Build the system
for i = 1:m+1
    w = c1 - c3/(1 + x(1,i)^2 + x(2,i)^2);
    x(1,i+1) = c4 + c2*(x(1,i)*cos(w) - x(2,i)*sin(w));
    x(2,i+1) = c2*(x(1,i)*sin(w) + x(2,i)*cos(w));
end
t_step = 1;
a = 0;
b = 200;
x = x(1,1:end-2);
N = size(x, 2);
y = awgn(x, noise_snr, 'measured');

fprintf('evaluating Gaussian kernel...')
tic
%Gaussian Kernel
KG = zeros(N);
kernel_size = median(pdist(transp(y))); %Gaussian kernel size sigma is median of pairwise distance of training data
p = 1/(2*kernel_size^2);
for i = 1:N
    for j = 1:N
        KG(i,j) = exp(-p*(y(i)-y(j)).^2);
    end
end

fprintf('evaluating KDS')

KV = zeros(N);
%t = (i-1)*t_step, tau = (j-1)*t_step
for i = 1:N
    for j = 1:N
        if j <= i
            KV(i,j) = a - (j)*t_step;
        else
            KV(i,j) = b - (j)*t_step;
        end
    end
end
KV = KV * 1/(b-a);

fprintf('evaluating B')

B = zeros(1,N);
C = zeros(N,N);
for i = 1:N
    B(i) = (-t_step / N) * sum(KG(:,i).*KV(:,i)) + (t_step^2 / N * 1/(b-a)) * sum(KG(:,i))*sum(KV(:,i));
end

fprintf('evaluating C')

for i = 1:N
    for j = 1:N
        %C(i,j) = sum(KV(:,i).*KV(:,j).*KG(i,j)) *t_step^2 / N;
        min_idx = min(i,j);
        max_idx = max(i,j);
        C(i,j) = (KV(1,i)*KV(1,j)*(min_idx-1) + KV(min_idx,i)*KV(min_idx,j)*(max_idx-min_idx) + KV(max_idx,i)*KV(max_idx,j)*(N-max_idx+1))*KG(i,j) *t_step^2 / N;
    end
end
fprintf('evaluating A')
A = -inv(C)*B';
%%
x_pred = zeros(size(y));
for i = 1:N
    for j = 1:N
        x_pred(i) = x_pred(i) + KDS_func(i, j, t_step, a, b, A(j))*y(j)*t_step;
    end
end

toc
error = x-x_pred;
%MSE(1, j) = sum((error(1,:).^2 + error(2,:).^2))/(m+1);
MSE = sum((error.^2))/N;
SNR = snr(x, error);

plot(y, 'Color', '#EDB120')
hold on
plot(x_pred, 'Color',	'#D95319')
plot(x , 'Color', '#0072BD')


legend('Noisy Signal', 'Estimated Signal','True Signal', 'Location','northwest')





