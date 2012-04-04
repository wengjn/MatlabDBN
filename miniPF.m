%This program is to implement particle filter with 3-state
%since it is a 3-state, you can not just use the one with example
%seems I got something wrong here

%clean the system first
clc
clear all
close all


%read data from mat file, get the true state var, without noise
dataFile = 'stateData.mat'; %contains true state variables
M = importdata(dataFile);
trueState = [];
for k=1:80
    trueState(1,k) = M(2,k);
    trueState(2,k) = M(3,k);
    trueState(3,k) = M(4,k);
end

%--------------------Prepare some variables------------------------------

Q = diag([0 0 0]); %System process noise covariance
R = 1; %measurement noise covariance
P = diag([0 0 0]); %initial estimation error covariance

x = [.1; .1; .1]; %initial state, x(1)=e3, x(2)=e8, x(3)=f10
xhat = [.1; .1; .2]; %initial state estimate

N = 200; %number of particles

%----------------------Initialize the particle fiter-----------------
for i = 1:N
    xhatplus(:, i) = x + sqrt(P) * [randn; randn; randn]; %standard particle filter
end

T = 0.5; %measurement timestep
randn('state', sum(100*clock)); %random number generator seed

tf = 80; %simulation length 
dt = 0.5; %timestep for integration
xArray = x;   
xhatArray = xhat; 

for t = T:T:tf
    %simulate the continuous system, with noise Q
    for tau = dt:dt:T
        %fourth order Runge Kutta integration
        [dx1, dx2, dx3, dx4] = RungeKutta(x);
        x = x + (dx1 + 2 * dx2 + 2 * dx3 + dx4)/6;
        x = x + sqrt(dt * Q) * [randn; randn; randn] * dt; %don't add Q
    end
    %simulate the noise measurement
    z = 1 + x(3) + sqrt(R) * randn;
    %simulate the continuous-time part of the particle filter (time update)
    xhatminus = xhatplus;
    for i = 1:N
        for tau = dt:dt:T
            %fourth order Runge Kutta integration
            %standard particle filter
            [dx1, dx2, dx3, dx4] = RungeKutta(xhatminus(:,i));
            xhatminus(:,i) = xhatminus(:,i) + (dx1 + 2 * dx2 + 2 * dx3 + dx4)/6;
            xhatminus(:,i) = xhatminus(:,i) + sqrt(dt * Q) * [randn; randn; randn] * dt;
        end
        zhat = 1 + xhatminus(3, i);
        vhat(i) = z - zhat;
    end
    %note that we need to scale all of the q(i) probabilities in a way that
    %does not change their relative magnitudes.
    %otherwise all of the q(i) elements will be zero because of the large
    %value of the exponential.
    %standard particle filter
    vhatscale = max(abs(vhat)) / 4;         
    qsum = 0;
    for i = 1:N
        q(i) = exp(-(vhat(i)/vhatscale)^2);    
        qsum = qsum + q(i);
    end
    
    %Normalize the likelihood of each a priori estimate
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
    xhatArray = [xhatArray xhat];
end


%draw the figures
close all

%t=1:tf;

figure;
plot(xArray(1,1:81));hold all
plot(xhatArray(1,1:81));
xlabel('timestep');
ylabel('state e4');
legend('True state', 'Particle filter estimate');

%plot(t, trueState(1, :));
% plot(t,trueState(1,:),'b.', t, xhatArray(1, :), 'k-');
% set(gca, 'FontSize',12); set(gcf, 'Color', 'White');
% xlabel('timestep');
% ylabel('state e4');
% legend('True state', 'Particle filter estimate');

figure;
plot(xArray(2,1:81));hold all
plot(xhatArray(2,1:81));
xlabel('timestep');
ylabel('state e8');
legend('True state', 'Particle filter estimate');
% plot(t,trueState(2,:), 'b.', t, xhatArray(2, :), 'k-');
% set(gca, 'FontSize', 12); set(gcf, 'Color', 'White');
% xlabel('timestep');
% ylabel('state e8');
% legend('True state', 'Particle filter estimate');

figure;
plot(xArray(3,1:81));hold all
plot(xhatArray(3,1:81));
xlabel('timestep');
ylabel('state f10');
legend('True state', 'Particle filter estimate');
% plot(t,trueState(3,:), 'b.', t, xhatArray(3, :), 'k-');
% set(gca, 'FontSize', 12); set(gcf, 'Color', 'White');
% xlabel('timestep');
% ylabel('state f10');
% legend('True state', 'Particle filter estimate');







% make them all in only one picture
% t = 1:80;
% figure;
% plot(t,e4Arr); hold all
% plot(t,e8Arr); hold all
% plot(t,f10Arr);
% xlabel('time step'); ylabel('state');
% legend('e4','e8','f10');
