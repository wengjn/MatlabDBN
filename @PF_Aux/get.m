function [x, P, pt, wt] = get(fltr)
% function [x, P, pt, wt] = get(fltr)
% input: 
% fltr: the filter object
% output: 
% x : the estimated mean
% P : the estimated covariance
% pt: the particles
% wt: the weights
%
% Author: Chihoon Lee, Lingji Chen
% Date  : Jan 20, 2006

wp = fltr.p .* repmat(fltr.w, size(fltr.p, 1), 1);
x = sum(wp, 2);
P = wp * fltr.p' - x * x';
pt = fltr.p;
wt = fltr.w;
