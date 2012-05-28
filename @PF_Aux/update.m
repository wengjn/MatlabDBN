function fltr = update(fltr, y)

% function fltr = update(fltr, y)
% input:
% fltr: the filter object
% y   : the current measurement
% output:
% fltr: the updated filter object
%
% Author: Lingji Chen
% Date  : March 15, 2006

%   Copyright (C) 2006, Lingji Chen, Chihoon Lee, Amarjit Budhiraja and 
%   Raman K. Mehra
%
%   This file is part of PFLib.
%
%   PFLib is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation; either version 2 of the License, or
%   (at your option) any later version.
%
%   PFLib is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with Foobar; if not, write to the Free Software
%   Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, 
%   MA  02110-1301  USA

% Bug fixes:
% June 11, 2007, by Charles Bartlett of The MITRE Corporation
%   changed the weight update line from 
%     fltr.w(i) = fltr.w(j) * density(fltr.v_d, y - y_pi)/ alpha(j);
%   to
%     fltr.w(i) = density(fltr.v_d, y - y_pi) / alpha(j);

N = length(fltr.w);
mu = fltr.p; 
alpha = zeros(1, N);

for i = 1:N
  mu(:, i) = fltr.f(fltr.p(:, i));
  yp = fltr.h(mu(:, i));
  alpha(i) = density(fltr.v_d, y - yp);
end;

wa = fltr.w .* alpha;
sumwa = sum(wa);
if sumwa <= realmin
  error('adjustment multiplier is numerically zero in Auxiliary PF');
end;
wa = wa / sumwa;
outIndex = fcn_ResampSimp(wa);

w_smpl = drawSamples(fltr.w_d, N);

for i = 1:N
  j = outIndex(i);
  fltr.p(:, i) = mu(:, j) + w_smpl(:, i);
  y_pi = fltr.h(fltr.p(:, i));
  fltr.w(i) = density(fltr.v_d, y - y_pi) / alpha(j);
end;

sum_w = sum(fltr.w);
if sum_w <= realmin
  error('weights are numerically zero; change parameters or method.');
end;
fltr.w = fltr.w / sum_w;


