%Auxiliary Particle Filter; 
close all; 
R = 1; % measurement noise variance (ft^2) 
Q = 10; % process noise covariance 
P = 1; % initial estimation error covariance 
 
x = .1; % initial state 
xhat = x; % initial state estimate 
 
N = 500; % number of particles   
 
% Initialize the particle filter. 
for i = 1 : N 
    xhatplusAux(i)=x + sqrt(P) * randn; 
    % auxiliary particle filter 
end 
 
tf = 50; % simulation length (seconds) 
 
xArray = x; 
yArray = [x^2/20+sqrt(P)*randn]; 
xhatAuxArray = xhat; 
 
for k = 2:tf 
        x = .5*x+25*x/(1+x^2)+8*cos(1.2*(k-1))+sqrt(Q)*randn; 
        y = x^2/20+sqrt(R)*randn; 
 
    for i = 1 : N 
        % auxiliary particle filter 
        xhatminusAux(i) = 0.5*xhatplusAux(i)+25*xhatplusAux(i)/(1+xhatplusAux(i)^2)+8*cos(1.2*(k-1))+sqrt(Q)*randn; 
        yhatAux = xhatminusAux(i)^2/20; 
        vhatAux(i) = y - yhatAux; 
    end 
    % Note that we need to scale all of the q(i) probabilities in a way 
    % that does not change their relative magnitudes. 
    % Otherwise all of the q(i) elements will be zero because of the 
    % large value of the exponential. 
     
      % auxiliary particle filter 
    vhatscaleAux = max(abs(vhatAux))/4 ; 
    qsumAux = 0; 
    for i = 1 : N 
        qAux(i) = exp(-(vhatAux(i)/vhatscaleAux)^2); 
        qsumAux = qsumAux + qAux(i); 
    end 
    % Regularize the probabilities - this is conceptually identical to the 
    % auxiliary particle filter - increase low probabilities and decrease 
    % high probabilities. 
    % Large k means low regularization (k = infinity is identical to the 
    % standard particle filter). Small k means high regularization (k = 1 
    % means that all probabilities are equal).  
    kAux = 100;  
    qAux = ((kAux - 1) * qAux + mean(qAux)) / kAux; 
    % Normalize the likelihood of each a priori estimate. 
    for i = 1 : N 
        qAux(i) = qAux(i) / qsumAux; 
    end 
    
    % Resample the auxiliary particle filter 
    for i = 1 : N 
        u = rand; % uniform random number between 0 and 1 
        qtempsum = 0; 
        for j = 1 : N 
            qtempsum = qtempsum + qAux(j); 
            if qtempsum >= u 
                xhatplusAux(i) = xhatminusAux(j); 
                break; 
            end 
        end 
    end 
    % The auxiliary particle filter estimate is the mean of the particles. 
    xhatAux = mean(xhatplusAux); 
    % Save data for plotting. 
    xArray = [xArray x]; 
    xhatAuxArray = [xhatAuxArray xhatAux]; 
end 
xhatAUX_RMS = sqrt((norm(xArray - xhatAuxArray))^2 / tf) 
figure(3);plot(xArray,'k-');hold on; plot(xhatAuxArray,'r-'); 
xlabel('time series'); ylabel('state');  
legend('True state', 'Auxiliary Particle filter estimate');  