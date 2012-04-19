%odetest test

clc
clear all

t0 = 5; %start
tf = 20; %stop
x0 = [1 -1 3]; %initial

[t, s] = ode23(@odetest, [t0, tf], x0);
x = s(:, 1);
y = s(:, 2);
z = s(:, 3);