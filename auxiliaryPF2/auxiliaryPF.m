%Implementation of Auxiliary Particle filter
%Refer to classic tutorial paper:
%Arulampalam. A tutorial on particle filters for online
%nonlinear/non-gaussian bayesian tracking. IEEE Transactions on Signal
%Processing. 50
%programming author: Jiannian Weng
%Date: May 15, 2012
%

function [xhAuxk, pfAux] = auxiliaryPF(yk, pfAux)
%note: right now I have bootstrap filter, genericPF and this auxiliaryPF
%see modifiedPF.m and particle_filter.m
%
%inputs:
%sysLC                  = system function handler
%observation            = observation function handler
%yk                     = observation vector at time step k
%pfAux                  = structure with the following fields
%     .k                = iteration number
%     .w                = weights (nx * Ns * T)
%     .Ns               = number of particles
%     .miuk             = characterization of x_k
%     .particles        = particles (nx * Ns * T)
%     .gen_x0           = function handler of a procedure that samples 
%                         initial pdf p_x0
%     .p_yk_given_xk    = function handler of the observation likelihood
%                         pdf p(y_k|x_k)
%     .gen_sys_noise    = function handler of a procedure that generate
%                         noise
%
%
%outputs:
%xhAuxk                 = estimated states
%pfAux                  = the same structure as in the input but updated
%                         at iteration k
%

k = pfAux.k;
if k == 1
    error('error: k must be an integer or equal than 2');
end

%initialize variables
Ns = pfAux.Ns; 
nx = size(pfAux.particles, 1);

wkm1 = pfAux.w(:,:, k-1);           %weights of laste iteration
if k==2
    for i = 1:Ns
        pfAux.particles(:, i, 1) = pfAux.gen_x0();
    end
    wkm1 = repmat(1/Ns, Ns, 1);
end

%memory space for particles
xkm1 = pfAux.particles(:, :, k-1);  %extract particles from last iteration
xk   = zeros(size(xkm1));
wk   = zeros(size(xkm1));
miuk = zeros(size(xkm1));

%iterations for Ns particles
%three step, first index k, then particle x, then weight for each particle
for i = 1:Ns
    %calculate miu_k_i first for every sample, sample from p(x_k | x_(k-1))
    miuk(i) = sysLC(k, xkm1(i), [pfAux.gen_sys_noise()]);
    
    wkm1(i);
    
    %sample index j from (71)
    %ji = wkm1(i) * pfAux.p_yk_given_miuk(k, yk, miuk(i));
    %ji = round(ji*100)

    while true
        ji = normrnd(0, 1);
        ji = round(abs(ji)*100);
        if ji >=1 && ji <= 100;
            break;
        end
    end
    
    %sample particles from q(x_k, i|y_k) using miu_k_i (68)
    %xk(:, i) = wkm1(i) * pfAux.p_yk_given_miuk(k, yk, miuk(:,i)) * miuk(:, i); % .* pfAux.p_xk_given_xkm1(k, xk(:,i), xkm1(:,i));
    xk(i) = sysLC(k, xkm1(ji), [pfAux.gen_sys_noise()]);
 
    
    %assign each sample with a weight ref to (72)
    wk(i) = pfAux.p_yk_given_xk(k, yk, xk(i)) ./ pfAux.p_yk_given_miuk(k, yk, miuk(i));
     % wk(:,i) = wkm1(i) * pfAux.p_yk_given_xk(k, yk, xk(:,i));
    % wk(:, i) = pfAux.p_yk_given_miuk(k, yk, miuk(:,i)) * wkm1(i);
    
end

%normalize weight vector
wk(1, :) = wk(1, :) ./ sum(wk(1, :));

%optional resampling or resample every time
Neff = 1 / sum(wk(1, :).^2);

resample_percent = 0.5;
Nt = resample_percent * Ns;
if Neff < Nt
    disp('Resampling ...');
    [xk, wk] = resample(xk, wk);
end


%compute estimated state
xhAuxk = zeros(nx, 1);
for i = 1:Ns
    xhAuxk = xhAuxk + wk(i) * xk(i);
end

%store new weights and particles
pfAux.w(:,:,k) = wk;
pfAux.particles(:,:,k) = xk;
pfAux.miuk(:,:,k) = miuk;

return;

%Resampling function
function [xk, wk] = resample(xk, wk)
Ns = length(wk);
c(1) = 0;
for i = 2:Ns
    c(i) = c(i-1) + wk(1, i);
end

for j = 1:Ns
    i = 1;
    u1 = (1/Ns) * rand();
    u(j) = u1 + (1/Ns) * (j-1);
    while u(j) >c(i)
        i = i +1 ;
        if i>=Ns
            break;
        end
        xk(1,j) = xk(1,i);
        wk(1,j) = 1/Ns;
    end
end

return;
    
