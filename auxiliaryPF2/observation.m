%observation equation to simulate the electrical LC system
%we have only one measurement

function [yk] = observation(k, xk, vk)

%linearize the states variable in this system
%vk is the observation noise


yk = xk^2/20 + vk;