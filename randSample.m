%Simulate the distribution observed in the 
%human random digit generation task

%example2: given a probabilities for each digits
%plot its pdf and cdf 
theta = [0.000;
         0.100;
         0.090;
         0.095;
         0.200;
         0.175;
         0.190;
         0.050;
         0.100;
         0.000; ];
     
seed = 1; rand('state', seed);

k = 1000;
digitset = 0:9;
y = randsample(digitset, k, true, theta);

figure(1); clf;

hist(y, digitset);
%counts = hist(y, digitset);
%bar(digitset, counts, 'k');
xlim([-0.5 9.5]);
xlabel('Digit');
ylabel('Frequency');
title('Distribution of simulated draws of human digit generator');

