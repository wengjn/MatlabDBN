%Implementation of Generic Particle filter
%Refer to classic tutorial paper:
%Arulampalam. A tutorial on particle filters for online
%nonlinear/non-gaussian bayesian tracking. IEEE Transactions on Signal
%Processing. 50
%programming author: Jiannian Weng
%Date: April 15, 2012
%

function [xhk, pf] = particle_filter(yk, pf)
%note: last time I wrote a script when resampling is performed on each step
%that is called the bootstrap particle filter, see modifiedPF.m
%
%inputs:
%sysLC       =  function handle to process equation
%observation =  function handle the measurement equation
%yk   =  observation vector at time k (column vector)
%pf   =  structure with the following fields
%  .k             = iteration number
%  .w             = weights (nx * Ns * T)
%  .Ns            = number of particles
%  .parcticles    = particles (nx * Ns * T)
%  .gen_x0        = function handle of a procedure that samples from the
%  initial pdf p_x0
%  .p_yk_given_xk = function handle of the observation likelihood pdf
%  p(y[k]|x[k])
%  .gen_sys_noise = function handle of a procedure that generates system
%  noise
%
%outputs:
%xhk = estimated state
%pf  = the same structure as in the input but updated at iteration k
%
%

k = pf.k;
if k==1
    error('error: k must be an iteger greater or equal than 2');
end

%initialize variables
Ns = pf.Ns;    %number of particles
nx = size(pf.particles, 1);   %number of states

wkm1 = pf.w(:,:, k-1);   %weights of last iteration
if k==2
    for i = 1:Ns
        pf.particles(:, i, 1) = pf.gen_x0(); %at time k=1
    end
    wkm1 = repmat(1/Ns, Ns, 1); %all particles have the same weight
end



%the importance sampling function:
%prior:  (I use this method here)
%q_xk_given_xkm1_yk = pf.p_xk_given_xkm1;

%optimal: (to be completed)
%q_xk_given_xkm1_yk = q_xk_given_xkm1^i_yk;
%note this pdf can be approximated by MCMC methods: they are expensive but
%they may be useful when non-iterative schemes fail

%memory space for particles
xkm1 = pf.particles(:,:,k-1); %extract particles from last iteration
xk   = zeros(size(xkm1));     % = zeros(nx, Ns);
wk   = zeros(size(xkm1));     % = zeros(nx, Ns);

%algorithm 3 of ref
for i = 1:Ns
    %xk(:,i) = sample_vector_from q_xk_given_xkm1_yk, given xkm1(:,i) and yk
    %using the prior pdf: pf.p_xk_given_xkm1: eq 62
    xk(:,i) = sysLC(k, xkm1(:,i), [0;0;pf.gen_sys_noise()]);
    
    %eq 48
    %wk(i) = wkm1(i) * p_yk_given_xk(yk,xk(:,i)) * p_xk_given_xkm1(xk(:,i),
    %xkm1(:,i))/q_xk_given_xkm1_yk(xk(:,i),xkm1(:,i),yk);
    
    %weights when using prior eq 63
    wk(:,i) = wkm1(i) * pf.p_yk_given_xk(k, yk, xk(:,i));
    
    %weights when using optimal eq 53
    %wk(i) = wkm1(i) * p_yk_given_xkm1(yk, xkm1(:,i)); %we don't know this
    %pdf
    
end

%normalize weight vector
%wk = wk ./ sum(wk);
wk(1,:) = wk(1,:) ./ sum(wk(1,:));
wk(2,:) = wk(2,:) ./ sum(wk(2,:));
wk(3,:) = wk(3,:) ./ sum(wk(3,:));


%calculate effective sample size: eq 51
Neff1 = 1/sum(wk(1,:).^2);

%Resampling
%remove this condition and sample on each iteration:
% [xk, wk] = resample(xk, wk, resampling_strategy);
%if you want to implement the bootstrap particle fileter
resample_percent = 0.50;
Nt = resample_percent*Ns;
if Neff1 < Nt
    disp('Resampling1...')
    [xk, wk] = resample(xk, wk);
    %{xk, wk} is an approximate discrete representation of p(x_k | y_{1:k}
end

Neff2 = 1/sum(wk(2,:).^2);
if Neff2 < Nt
    disp('Resampling2...')
    [xk, wk] = resample2(xk, wk);
end

Neff3 = 1/sum(wk(3,:).^2);
if Neff3 < Nt
    disp('Resampling3...')
    [xk, wk] = resample3(xk, wk);
end

%compute estimated state
xhk = zeros(nx, 1);
for i = 1:Ns
    xhk = xhk + wk(i)*xk(:,i);
end

%store new weights and particles
pf.w(:,:,k) = wk;
pf.particles(:,:,k) = xk;

return;

%REsampling function
function [xk, wk] = resample(xk, wk)
Ns = length(wk);
c(1) = 0; %initialize cdf
for i = 2:Ns
    %construct cdf
    c(i) = c(i-1) + wk(1,i);
end


for j = 1:Ns
    %start from bottom of cdf
    i = 1;
    %draw a starting point
    u1 = (1/Ns)*rand();
    u(j) = u1 + (1/Ns) * (j - 1);
    while u(j) > c(i)
        i = i + 1;
        if i >= Ns
            break;
        end
        
    end
    xk(1, j) = xk(1, i);
    wk(1, j) = 1/Ns;
end

return;

function [xk, wk] = resample2(xk, wk)
Ns = length(wk);
c(1) = 0; %initialize cdf
for i = 2:Ns
    %construct cdf
    c(i) = c(i-1) + wk(2,i);
end


for j = 1:Ns
    %start from bottom of cdf
    i = 1;
    %draw a starting point
    u1 = (1/Ns)*rand();
    u(j) = u1 + (1/Ns) * (j - 1);
    while u(j) > c(i)
        i = i + 1;
        if i >= Ns
            break;
        end
        
    end
    xk(2, j) = xk(2, i);
    wk(2, j) = 1/Ns;
end

return;


function [xk, wk] = resample3(xk, wk)
Ns = length(wk);
c(1) = 0; %initialize cdf
for i = 2:Ns
    %construct cdf
    c(i) = c(i-1) + wk(3,i);
end


for j = 1:Ns
    %start from bottom of cdf
    i = 1;
    %draw a starting point
    u1 = (1/Ns)*rand();
    u(j) = u1 + (1/Ns) * (j - 1);
    while u(j) > c(i)
        i = i + 1;
        if i >= Ns
            break;
        end
        
    end
    xk(3, j) = xk(3, i);
    wk(3, j) = 1/Ns;
end

return;














