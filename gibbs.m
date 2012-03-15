%This is an example of Gibbs Sampling

N = 10000;
n = 10;
c = 15;
X = 2*ones(1, n);

for m = 1:N; i = ceil(n*rand);
    S = sum(X) - X(i);
    X(i) = max(c-S, 0) - log(rand)/i;
    H(m) = S + X(i);
end

disp([mean(H) 2*std(H)/sqrt(N)]);

N = 10000000; XR = -diag(1./[1:n])*log(rand(n,N));
S = sum(XR); H = S.*(S>15)/mean(S>15);
disp([mean(H) 2*std(H)/sqrt(N)]);