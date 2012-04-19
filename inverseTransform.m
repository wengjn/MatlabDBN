%inversion transform sampling algorithm
%The code here is to compute the integration of y = xexp(x)dx [0 1]

N = 100;
for i = 1:N
    u = rand;
    x(i) = log(1 + (exp(1)-1)*u);
end

y = (exp(1) - 1)*sum(x) / N;
