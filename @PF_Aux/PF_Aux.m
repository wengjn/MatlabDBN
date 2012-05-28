function fltr = PF_Aux(sys, init_distr, opts)
% function fltr = PF_Aux(sys, init_distr, opts)
% PF_Aux: class constructor for an Auxiliary  Particle Filter, with a simple
% proposal (the state transition) but modified resampling.
% 
% sys       : a SigObsSys object
% init_distr: initial distribution of the particles
% opts      : a structure array of options, obtained from calling setPFOptions
% fltr      : the contructed filter
%
% Author: Lingji Chen
% Date: March 15, 2006


if nargin == 0
  fltr.f = []; fltr.h = []; fltr.w_d = []; fltr.v_d = [];
  fltr.N = []; 
  fltr.p = []; fltr.w = []; 
  fltr = class(fltr,'PF_Aux');
elseif isa(sys, 'PF_Aux')
  fltr = sys;
elseif nargin ~= 3
  error('wrong number of arguments');
else
  if ~isa(sys, 'SigObsSys')
    error('wrong first argument: must be a SigObsSys object');
  end;
  if ~setPFOptions(opts)
    error(['wrong second argument: must be a valid option' ...
	   ' structure']);
  end;
  [x, y, f, h, w_d, v_d] = get(sys);
  fltr.f = f; % state equation
  fltr.h = h; % measurement equation
  fltr.w_d = w_d; % distribution of state noise
  fltr.v_d = v_d; % distribution of measurement noise
  try
    fltr.N = opts.NumPart;
  catch
    error(['wrong option structuure:' lasterr]);
  end;
  
  try
    fltr.p = drawSamples(init_distr, fltr.N); % particles
    fltr.w = ones(1, fltr.N) / fltr.N;        % weights
  catch
    error('wrong initial particle distribution');
  end;
  fltr = class(fltr, 'PF_Aux');
end
