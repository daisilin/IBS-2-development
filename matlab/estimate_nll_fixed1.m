function [nll,nll_sd,output]=estimate_nll_fixed1(model,stim,resp_real,theta,Nsamples,Nreps)
%ESTIMATE_NLL_FIXED1 Negative log likelihood estimation via fixed sampling.
% (Uses correction of one pseudo-count. prior added: assume one hit and one miss have been seen;)

persistent samples_used;
persistent reps_used;
persistent funcalls;
persistent p_vec;

if isempty(samples_used)
    samples_used = 0;
    reps_used = 0;
    funcalls = 0;
    p_vec = 0;
end
if nargin < 4 || isempty(theta)
    nll = []; nll_sd = []; 
    output.samples_used = samples_used;
    output.reps_used = reps_used;
    output.funcalls = funcalls;
    output.p_vec = p_vec; 
    samples_used = 0;
    reps_used = 0;
    funcalls = 0;
    p_vec = 0;
    return;
end

if nargin < 6 || isempty(Nreps); Nreps = 1; end

Ntrials = size(resp_real,1);
m = zeros(Ntrials,1);

Neff = Nsamples*Nreps; %number of fixed samples * number of repeats
for iSample = 1:Neff
    resp = feval(['generate_resp_' model],stim,theta,1:Ntrials);
    m = m + all(resp==resp_real,2); % m + a vector oflogical where hits = 1;m counts the number of hits for each trial
end

% Estimate p with one pseudo-count correction, a vector
p = (m+1)/(Neff+2); % add prior: assume one hit and one miss have been seen; denomenator changed from +1 to +2

nll = -sum(log(p));

if nargout > 1
    % Estimate SD via bootstrap
    p_eff = m/Neff;    
    nll_bootstrap = [];    
    Niters=5;
    Nbootperiter=10;  
    for iBoot = 1:Niters
        m_bootstrap = binornd(Neff,repmat(p_eff,[1,Nbootperiter]));
        p_bootstrap = (m_bootstrap+1)/(Neff+1);
        nll_bootstrap = [nll_bootstrap, -sum(log(p_bootstrap),1)];
    end
    nll_sd = std(nll_bootstrap);
end

samples_used = samples_used + Neff;
funcalls = funcalls + 1;
reps_used = reps_used + Nreps;
p_vec = p;  %add persistent variable, a vector of likelihood 

if nargout > 2
    output.samples_used = samples_used;
    output.reps_used = reps_used;    
    output.funcalls = funcalls;
    output.p_vec = p_vec;
end
