%visualize normal distribution

%example1: explore the normal distribution N(mu, sigma)
mu = 100;
sigma = 15;
xmin = 70;
xmax = 130;
n = 100;
k = 10000;

x = linspace(xmin, xmax, n);
p = normpdf(x, mu, sigma);
c = normcdf(x, mu, sigma);

figure(1); clf;

subplot(1, 3, 1);
plot(x, p, 'k-');
xlabel('x'); ylabel('pdf');
title('Probability Density Function');

subplot(1,3,2);
plot(x, c, 'k-');
xlabel('x'); ylabel('cdf');
title('Cumulative Density Function');

y = normrnd(mu, sigma, k, 1);
subplot(1,3,3);
hist(y, 20);
xlabel('x'); ylabel('frequency');
title('Histogram of random values');



