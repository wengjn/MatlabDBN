%MC method applications

%Transformation of pdf's
%step1: sample from uniform distribution, theta
%step2: get the invert of F(x)
%step3, get the sample F^(-1)(theta)
N = 2000
xArr = [];
for t = 1:N
    theta = rand;
    x = -.5 * log(1-theta);
    xArr = [xArr x];
end
tf = 1:1:N;
%plot(tf, xArr, 'b-')

%I=Pr(Z>3) when Z~N(0,1)
%Monte Carlo method to estimate a probability
%The idea is to sample from normal distribution and to count the number of 
%observations that is greater than 3
format long
for i = 1:100
    a(i) = sum(randn(100,1)>=3)/100;
end
meanMC = mean(a);
varMC = var(a);

%using importance sampling method
%http://www.wikicoursenote.com/wiki/Importance_Sampling_and_Markov_Chain_Monte_Carlo_(MCMC)

for j = 1:100
    N = 100;
    x = randn(N, 1) + 4;
    for ii = 1:N
        h = x(ii)>=3;
        b = exp(8-4*x(ii));
        w(ii) = h*b;
    end
    I(j) = sum(w)/N;
end
MEAN = mean(I);
VAR = var(I);

%
x = -200:200;
for k = 1:401
    y(k) = (1/sqrt(2*pi))*exp(-x(k)*x(k)/2);
    m(k) = (1/sqrt(2*pi)/4)*exp(-x(k)^2/2/16);
end
plot(x,y)
hold all
plot(x,m)