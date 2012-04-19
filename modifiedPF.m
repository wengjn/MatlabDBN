%This program is to implement particle filter with 3-state
%since it is a 3-state, you can not just use the one with example
%rightLC.mdl simulate such system in simulink
%you have to use variable-step size ode45 to solve such system,
%otherwise, you can not see that is is oscilating, too small

%clean the system first
clc
clear all
close all


%read data from mat file, get the true state var, without noise
dataFile = 'stateData.mat'; %contains true state variables
M = importdata(dataFile);
trueState = [];
for k=1:90
    trueState(1,k) = M(2,k);
    trueState(2,k) = M(3,k);
    trueState(3,k) = M(4,k);
end

%--------------------Prepare some variables------------------------------

Q = diag([0 0 1.0e-6]); %System process noise covariance1.0e-6
R = 2; %measurement noise covariance
P = diag([0 0 0]); %initial estimation error covariance

x = [0; 0; 0]; %initial state, x(1)=e4, x(2)=e8, x(3)=f10
xhat = [0.5; 0; .2]; %initial state estimate

N = 200; %number of particles

%----------------------Initialize the particle fiter-----------------
for i = 1:N
    xhatplus(:, i) = x + sqrt(P) * [randn; randn; randn]; %standard particle filter
end

h = 0.15; %measurement timestep
%randn('state', sum(100*clock)); %random number generator seed

tf = 500; %simulation length 
dt = 0.5; %timestep for integration
xArray = x;  %with noise state
xtArray = x; %without noise, true state
xhatArray = xhat; 

for t = 0:h:tf
    %simulate the continuous system, with noise Q
    %fourth order Runge Kutta integration
    [slope] = EulerMethod(x);
    x = x + h*slope;
    xt = x;
    x = x + sqrt(dt * Q) * [randn; randn; randn] * dt; %make Q=0 at first
    
    %simulate the noise measurement
    %z = 1 + x(3) + sqrt(R) * randn;
    z = x(1) + 2*x(2) + 3*x(3) + sqrt(R)*randn;
    %simulate the continuous-time part of the particle filter (time update)
    xhatminus = xhatplus;
    for i = 1:N
            %fourth order Runge Kutta integration
            %standard particle filter
        [slope] = EulerMethod(xhatminus(:,i));
        xhatminus(:,i) = xhatminus(:,i) + h*slope;
        xhatminus(:,i) = xhatminus(:,i) + sqrt(dt * Q) * [randn; randn; randn] * dt;
        %zhat = 1 + xhatminus(3, i);
        zhat = xhatminus(1,i) + 2*xhatminus(2, i) + 3*xhatminus(3, i);
        vhat(i) = z - zhat;
        q(i) = (1/sqrt(R)/sqrt(2*pi))*exp(-vhat(i)^2/2/R); %N~(0, R)
    end
    %note that we need to scale all of the q(i) probabilities in a way that
    %does not change their relative magnitudes.
    %otherwise all of the q(i) elements will be zero because of the large
    %value of the exponential.
    %standard particle filter
    
    %Normalize the likelihood of each a priori estimate
    qsum = sum(q); %sum all the weights
    for i = 1:N
        q(i) = q(i) / qsum;
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
    %save data for plotting
    xArray = [xArray x];
    xtArray = [xtArray xt];
    xhatArray = [xhatArray xhat];
end


%draw the figures
close all

t = 0:h:tf;

figure;
plot(t, xArray(1,1:length(t)),'b-');figure(gcf)
hold all 
plot(t, xhatArray(1,1:length(t)),'k-');figure(gcf)
hold all
plot(t, xtArray(1,1:length(t)),'r-');
axis normal
xlabel('timestep');
ylabel('state e4');
legend('True state with noise', 'Particle filter estimate', 'without noise');


figure;
plot(t, xArray(2,1:length(t)),'b-');figure(gcf)
hold all 
plot(t, xhatArray(2,1:length(t)),'k-');figure(gcf)
hold all
plot(t, xtArray(2,1:length(t)),'r-');
axis normal
xlabel('timestep');
ylabel('state e8');
legend('True state with noise', 'Particle filter estimate', 'without noise');

figure;
plot(t, xArray(3,1:length(t)),'b-');figure(gcf)
hold all 
plot(t, xhatArray(3,1:length(t)),'k-');figure(gcf)
hold all
plot(t, xtArray(3,1:length(t)),'r-');
axis normal
xlabel('timestep');
ylabel('state f10');
legend('True state with noise', 'Particle filter estimate','without noise');



