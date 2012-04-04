%kalman filter

clear all;
close all;

%initial parameters
n_iter = 50;
sz = [n_iter, 1];
x = -0.37727;
z = x + sqrt(0.01)*randn(sz);

Q = 1e-5;

%allocate space for arrays
xhat = zeros(sz);
P = zeros(sz);
xhatminus = zeros(sz);
Pminus = zeros(sz);
K = zeros(sz);

R = 0.01;

%initial guesses
xhat(1) = 0.0;
P(1) = 1.0;

for k = 2:n_iter
    %time update
    xhatminus(k) = xhat(k-1);
    Pminus(k) = P(k-1) + Q;
    %measument update
    K(k) = Pminus(k)/(Pminus(k)+R);
    xhat(k) = xhatminus(k) + ...
        K(k)*(z(k)-...
        xhatminus(k));
    P(k) = (1-K(k))*Pminus(k);
end
figure();
plot(z, 'k+');
hold on;
plot(xhat, 'b-');
hold on;
plot([0:0.1:50],x,'g-');