%This is a case study using exactly as the classic tutorial paper 
%to solve electrical system
%
%programming author: Jiannian Weng
%Date: April 15, 2012
%

clear all
close all
clc

%process equation
nx = 3;                         %number of states, voltage on c1, c2 and current through L
%see the function sysLC()

%observation equation
%see the function observation()
ny = 1;                         %number of observations, only one measurement

%PDF of process noise and noise generator function
nu = 3;                         %size of the vector of process noise
sigma_u = sqrt(1.0e-6);         %Q
p_sys_noise = @(u) normpdf(u, 0, sigma_u);
gen_sys_noise = @(u) normrnd(0, sigma_u);    %sample from noise


%PDF of observation noise and noise generator function
nv = 1;                         %size of the vector of observation noise
sigma_v = sqrt(2);              %R
p_obs_noise = @(v) normpdf(v, 0, sigma_v);
gen_obs_noise = @(v) normrnd(0, sigma_v);    %sample from noise


%Initial PDF
p_x0 = @(x) normpdf(x, 0, sqrt(2));
gen_x0 = @(x) normrnd(0, sqrt(2));           %sample from p_x0

%Transition prior PDF p(x[k]|x[k-1])
p_xk_given_xkm1 = @(k, xk, xkm1) p_sys_noise(xk - sysLC(k, xkm1, 0));

%Observation likelihood PDF p(y[k]|x[k])
p_yk_given_xk = @(k, yk, xk) p_obs_noise(yk - observation(k, xk, 0));

%number of time steps
T = 500;

%memory space for true state
x = zeros(nx, T); y = zeros(ny, T);
u = zeros(nu, T); v = zeros(nv, T);

%simulate system
xh0 = [0; 0; 0];                            %column vector
u(:, 1) = [0;0;gen_sys_noise()];
v(:, 1) = gen_obs_noise(sigma_v);
x(:, 1) = xh0;
y(:, 1) = observation(1, xh0, v(:,1));

for k = 2:T
    %here we are basically sampling from p_xk_given xkm1 and from
    %p_yk_given_xk
    u(:, k) = [0;0;gen_sys_noise()];
    v(:, k) = gen_obs_noise();
    x(:, k) = sysLC(k, x(:, k-1), u(:, k));     %simulate state
    y(:, k) = observation(k, x(:, k-1), v(:, k)); %simulate measurement
end


%memory space for particle filter estimate
xh = zeros(nx, T); xh(:, 1) = xh0;
yh = zeros(ny, T); yh(:, 1) = observation(1, xh0, 0);

pf.k = 1;                           %initial iteration
pf.Ns = 100;                        %number of particles
pf.w = zeros(nx, pf.Ns, T);         %weights
pf.particles = zeros(nx, pf.Ns, T); %particles
pf.gen_x0 = gen_x0;                 %function for sampling from initial pdf
pf.p_yk_given_xk = p_yk_given_xk;   %function of the observation likelihood pdf
pf.gen_sys_noise = gen_sys_noise;   %function for generating system noise

%estimate state using particle filter
for k = 2:T
    fprintf('Iteration = %d/%d\n', k, T);
    %state estimation
    pf.k = k;
    [xh(:, k), pf] = particle_filter(y(:, k), pf);
    %filtered observation
    yh(:, k) = observation(k, xh(:, k), 0);
end


%plots
figure 
hold on;
h1 = plot(1:T, y, 'r');
h2 = plot(1:T, yh, 'b');
legend([h1 h2], 'observation', 'filtered observation');
title('Observation vs filtered observation by the particle filter');

figure 
hold on;
h1 = plot(1:T, x(3,:), 'r');
h2 = plot(1:T, xh(3, :), 'b');
legend([h1,h2], 'true state3', 'filtered state3');
title('true state vs filtered state3 by particle filter');
return;







