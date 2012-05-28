function [StdRMSErr, AuxRMSErr] = ParticleFilter

clear all;
close all;
clc;

%compare the generic particle filter with auxiliary particle filter


Q = diag([0 0 1.0e-6]); %System process noise covariance1.0e-6
R = 2; %measurement noise covariance
P = diag([0 0 0]); %initial estimation error covariance

x = [0; 0; 0]; %initial state, x(1)=e4, x(2)=e8, x(3)=f10
xhat = [0.5; 0; .2]; %initial state estimate

N = 200; %number of particles


%----------------------Initialize the particle fiter-----------------
for i = 1:N
    xhatplus(:, i) = x + sqrt(P) * [randn; randn; randn]; %generic particle filter
    xhatplusAux(:, i) = xhatplus(:, i); %auxiliary particle filter
end

h = 0.15; %measurement timestep

tf = 200; %simulation length 
dt = 0.5; %timestep for integration
xArray = x;  %with noise state
xtArray = x; %without noise, true state
xhatArray = xhat; 
xhatAuxArray = xhat;
yhArray =  [0];
yhAuxArray = [0];

for t = 0:h:tf
    %simulate the continuous system, with noise Q
    %fourth order Runge Kutta integration
    fprintf('Particle Filter Iteration = %d\n', t);
    [slope] = EulerMethod(x);
    x = x + h*slope;
    xt = x;
    x = x + sqrt(Q) * [randn; randn; randn]; %make Q=0 at first
    
    %simulate the noise measurement
    z = x(1) + 2*x(2) + 3*x(3) + sqrt(R)*randn;
    %simulate the continuous-time part of the particle filter (time update)
    xhatminus = xhatplus;
    xhatminusAux = xhatplusAux;
    for i = 1:N
            %standard particle filter
        [slope] = EulerMethod(xhatminus(:,i));
        xhatminus(:,i) = xhatminus(:,i) + h*slope;
        xhatminus(:,i) = xhatminus(:,i) + sqrt( Q) * [randn; randn; randn];
        %zhat = 1 + xhatminus(3, i);
        zhat = xhatminus(1,i) + 2*xhatminus(2, i) + 3*xhatminus(3, i);
        vhat(i) = z - zhat;
        q(i) = (1/sqrt(R)/sqrt(2*pi))*exp(-vhat(i)^2/2/R); %N~(0, R)
        
            %auxiliary particle filter
        [slope] = EulerMethod(xhatminusAux(:, i));
        xhatminusAux(:, i) = xhatminusAux(:, i) + h*slope;
        xhatminusAux(:, i) = xhatminusAux(:, i) + sqrt( Q) * [randn; randn; randn];
        zhatAux = xhatminusAux(1, i) + 2*xhatminusAux(2, i) + 3*xhatminusAux(3, i);
        vhatAux(i) = z - zhatAux;
    end
    %note that we need to scale all of the q(i) probabilities in a way that
    %does not change their relative magnitudes.
    %otherwise all of the q(i) elements will be zero because of the large
    %value of the exponential.
    %standard particle filter
    
    %Normalize the likelihood of each a priori estimate
    vhatscale = max(abs(vhat)) / 4;
    qsum = 0;
    for i = 1:N
        q(i) = exp(-(vhat(i)/vhatscale)^2);
        qsum = qsum + q(i);
    end
    for i = 1:N
        q(i) = q(i)/qsum;
    end
    
    %auxiliary particle filter
    vhatscaleAux = max(abs(vhatAux)) / 4;
    qsumAux = 0;
    for i = 1:N
        qAux(i) = exp(-(vhatAux(i)/vhatscaleAux)^2);
        qsumAux = qsumAux + qAux(i);
    end
    kAux = 1.1;
    qAux = ((kAux - 1) * qAux + mean(qAux)) / kAux;
    for i = 1:N
        qAux(i) = qAux(i) / qsumAux;
    end
    
    %Resample the standard particle filter
    for i = 1:N
        u = rand; %uniform random number between 0 and 1
        qtempsum = 0;
        for j = 1:N
            qtempsum = qtempsum + q(j);
            if qtempsum >= u
                xhatplus(:,i) = xhatminus(:,j);
                break;
            end 
        end
    end
    %The standard particle filter estimate is the mean of the particles.
    xhat = mean(xhatplus')';
    
    %Resample the auxiliary particle filter
    for i = 1:N
        u = rand;
        qtempsum = 0;
        for j = 1:N
            qtempsum = qtempsum + qAux(j);
            if qtempsum >= u
                xhatplusAux(:,i) = xhatminusAux(:,j);
                break;
            end
        end
    end
    xhatAux = mean(xhatplusAux')';
    
    yh = xhat(1) + 2*xhat(2) + 3*xhat(3);
    yhAux = xhatAux(1) + 2*xhatAux(2) + 3*xhatAux(3);
    
    %save data for plotting
    xArray = [xArray x];
    xtArray = [xtArray xt];
    xhatArray = [xhatArray xhat];
    xhatAuxArray = [xhatAuxArray xhatAux];
    yhArray = [yhArray yh];
    yhAuxArray = [yhAuxArray yhAux];
end

%draw the figures
close all

t = 0:h:tf;
figure;
semilogy(t, abs(xArray(1, 1:length(t)) - xhatArray(1, 1:length(t))), 'b-'); hold;
semilogy(t, abs(xArray(1, 1:length(t)) - xhatAuxArray(1, 1:length(t))), 'r:');
set(gca, 'FontSize', 12); set(gcf, 'Color', 'White');
xlabel('Seconds');
ylabel('state e4 estimation error');
legend('Generic particle filter', 'Auxiliary particle filter');

figure;
plot(t, yhArray(1:length(t)), 'b-'); figure(gcf); hold;
plot(t, yhAuxArray(1:length(t)), 'r-'); figure(gcf); hold;
xlabel('time step');
ylabel('observation');
legend('GPF observ', 'APF observ');

% 
% figure;
% plot(t, xArray(1,1:length(t)),'b-');figure(gcf)
% hold all 
% plot(t, xhatArray(1,1:length(t)),'k-');figure(gcf)
% hold all
% plot(t, xtArray(1,1:length(t)),'r-');
% hold all
% plot(t, xhatAuxArray(1,1:length(t)), 'g-');figure(gcf)
% axis normal
% xlabel('timestep');
% ylabel('state e4');
% legend('True state with noise', 'Particle filter estimate', 'without noise', 'Auxiliary PF');
% 
% 
% figure;
% plot(t, xArray(2,1:length(t)),'b-');figure(gcf)
% hold all 
% plot(t, xhatArray(2,1:length(t)),'k-');figure(gcf)
% hold all
% plot(t, xtArray(2,1:length(t)),'r-');
% axis normal
% xlabel('timestep');
% ylabel('state e8');
% legend('True state with noise', 'Particle filter estimate', 'without noise');
% 
% figure;
% plot(t, xArray(3,1:length(t)),'b-');figure(gcf)
% hold all 
% plot(t, xhatArray(3,1:length(t)),'k-');figure(gcf)
% hold all
% plot(t, xtArray(3,1:length(t)),'r-');
% axis normal
% xlabel('timestep');
% ylabel('state f10');
% legend('True state with noise', 'Particle filter estimate','without noise');
% 



