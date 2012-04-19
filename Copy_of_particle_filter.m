%Generic Particle filter

function [xhk, pf] = particle_filter(sys, yk, pf, resampleing_strategy)
%note: when resampling is performed on each step this algorithm is called
% the bootstrap particle filter

%inputs:
%sys  =  function handle to process equation
%yk   =  observation vector at time k (column vector)
%pf   =  structure with the following fields
%  .k             = iteration number
%  .w             = weights (Ns * T)
%  .Ns            = number of particles
%  .parcticles    = particles (nx * Ns * T)
%  .gen_x0        = function handle of a procedure that samples from the
%  initial pdf p_x0
%  .p_yk_given_xk = function handle of the observation likelihood pdf
%  p(y[k]|x[k])
%  .gen_sys_noise = function handle of a procedure that generates system
%  noise
%resampling_strategy = resampling strategy. set it either to
%'multinomial_resampling' or 'systematic_resampling'
%
%outputs:
%xhk = estimated state
%pf  = the same structure as in the input but updated at iteration k
%
%reference:
%Arulampalam. A tutorial on particle filters for online
%nonlinear/non-gaussian bayesian tracking. IEEE Transactions on Signal
%Processing. 50
%
%programming author: Jiannian Weng
%April 17, 2012

%
k = pf.k;
if k==1
    error('error: k must be an iteger greater or equal than 2');
end

%initialize variables
Ns = pf.Ns;    %number of particles
nx = size(pf.particles, 1);   %number of states

wkm1 = pf.w(:, k-1);   %weights of last iteration
if k==2
    for i = 1:Ns
        pf.particles(:, i, 1) = pf.gen_x0(); %at time k=1
    end
    wkm1 = repmat(1/Ns, Ns, 1);  %all particles have the same weight
end



%the importance sampling function:
%prior: (this method is sensitive to outliers)
%q_xk_given_xkm1_yk = pf.p_xk_given_xkm1;

%optimal:
%q_xk_given_xkm1_yk = q_xk_given_xkm1^i_yk;
%note this pdf can be approximated by MCMC methods: they are expensive but
%they may be useful when non-iterative schemes fail

%separate memory
xkm1 = pf.particles(:,:,k-1); %extract particles from last iteration
xk   = zeros(size(xkm1));     % = zeros(nx, Ns);
wk   = zeros(size(xkm1));     % = zeros(Ns, 1);

%algorithm 3 of ref
for i = 1:Ns
    %xk(:,i) = sample_vector_from q_xk_given_xkm1_yk, given xkm1(:,i) and yk
    %using the prior pdf: pf.p_xk_given_xkm1: eq 62
    xk(:,i) = sys(k, xkm1(:,i), pf.gen_sys_noise());
    
    %eq 48
    %wk(i) = wkm1(i) * p_yk_given_xk(yk,xk(:,i)) * p_xk_given_xkm1(xk(:,i),
    %xkm1(:,i))/q_xk_given_xkm1_yk(xk(:,i),xkm1(:,i),yk);
    
    %weights when using prior eq 63
    wk(i) = wkm1(i) * pf.p_yk_given_xk(k,yk,xk(:,i));
    
    %weights when using optimal eq 53
    %wk(i) = wkm1(i) * p_yk_given_xkm1(yk, xkm1(:,i)); %we don't know this
    %pdf
    
end

%normalize weight vector
wk = wk./sum(wk);

%calculate effective sample size: eq 51
Neff = 1/sum(wk.^2);

%Resampling
%remove this condition and sample on each iteration:
% [xk, wk] = resample(xk, wk, resampling_strategy);
%if you want to implement the bootstrap particle fileter
resample_percentaje = 0.50;
Nt = resample_percentaje*Ns;
if Neff < Nt
    disp('Resampling...')
    [xk, wk] = resample(xk, wk, resampleing_strategy);
    %{xk, wk} is an approximate discrete representation of p(x_k | y_{1:k}
end

%compute estimated state
xhk = zeros(nx, 1);
for i = 1:Ns
    xhk = xhk + wk(i)*xk(:,i);
end

%store new weights and particles
pf.w(:,k) = wk;
pf.particles(:,:,k) = xk;

return;

%REsampling function
function [xk, wk, idx] = resample(xk, wk, resampleing_strategy)

Ns = length(wk);  %Ns = number of particles
%wk = wk./sum(wk); %normalize weight vector, already done
switch resampleing_strategy
    case 'multinomial_resampling'
        with_replacement = true;
        idx = randsample(1:Ns, Ns, with_replacement, wk);
    case 'systematic_resampling'
        edges = min([0 cumsum(wk)]', 1); %protect against accumulated roundoff
        edges(end) = 1; %get the upper edge exact
        u1 = rand/Ns;
        [~, idx] = histc(u1:1/Ns:1, edges);
    otherwise
        error('Resampling strategy not implemented')
end

xk = xk(:,idx);             %extract new particles
wk = repmat(1/Ns, 1, Ns);   %now all particles have the same weight

return 

















