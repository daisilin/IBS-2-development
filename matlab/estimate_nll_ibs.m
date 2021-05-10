function [nll,nll_sd,output]=estimate_nll_ibs(model,stim,resp_real,theta,Nreps,thresh)
%ESTIMATE_NLL_IBS Negative log likelihood estimation via inverse binomial sampling.
% %thresh is for early stopping, if nll of data for each ibs run/repeat is smaller than a number; 
persistent samples_used;
persistent reps_used;
persistent funcalls;
%persistent K_tot; % trials x funcalls; nll for each trial, with all repeats averaged for each trial; each column is a different x (param setting)depending on bads
persistent allocate_reps; % num of reps that should be allocated for this x

if isempty(samples_used)
    samples_used = 0;
    reps_used = 0;
    funcalls = 0;
    %K_tot = [];
    allocate_reps = 0;
end
if nargin < 4 || isempty(theta)
    nll = []; nll_sd = []; 
    output.samples_used = samples_used;     
    output.reps_used = reps_used;
    output.funcalls = funcalls;
    %output.K_tot = K_tot;
    output.allocate_reps = allocate_reps; 
    samples_used = 0;
    reps_used = 0;
    funcalls = 0;
    return;
end

if nargin < 6 || isempty(thresh); thresh = Inf; end

Nreps_max = max(Nreps);

nll_vec = zeros(1,Nreps_max); %use the largest repeat num; Nreps is a vector of length = N repeats,
if nargout > 1; nll_var_vec = zeros(1,Nreps_max); end

Ntrials = size(resp_real,1);
if isscalar(Nreps); Nreps = Ntrials*ones(Ntrials,1); end

for iRep = 1:Nreps_max %loop over repeats
    resp = NaN(Ntrials,1);
    tries = zeros(Ntrials,1); % vector that counts how many tries/miss there have been before a match/hit for all trials; it gets reset for every repeat
    ind = iRep<=Nreps;
    n = sum(ind); %number of "active" trials (haven't matched) left
    tries(~ind) = NaN;
    nll = 0;
    
    % Note: thresh here only makes sense in iterations that use all trials
    while n>0 && nll<thresh % while there are trials that haven't been matched && nll of data for this repeat is smaller than a threshold
        resp(ind) = feval(['generate_resp_' model],stim,theta,ind); %generate synthetic response from the model; subsequent resps are for miss (and repeated) trials only
        ind = any(resp~=resp_real,2); %find trial index where true and synthetic responses are not matched; logical array where unmatched trials are set to 1
        tries = tries+ind; %by adding ind, a logical vector where unmatched tirals = 1 ; we get a vector with Ntiral length, each idx counts how many misses for each (repeated) trial
        n = sum(ind);
        
        K = tries+1;
        Ktab = psi(1:max(K(:)))' - psi(1);        
        nll = nansum(Ktab(K));
    end
    
    nll_vec(iRep) = nansum(Ktab(K)./Nreps);    
    
    %Compute estimate of the variance if requested
    if nargout > 1
        Ktab = -(psi(1,1:max(K(:)))' - psi(1,1));
        nll_var_vec(iRep) = nansum(Ktab(K)./Nreps.^2);
    end
    
    samples_used = samples_used + nanmean(tries+1);
end

funcalls = funcalls + 1;
reps_used = reps_used + mean(Nreps);

nll = sum(nll_vec); % return the average nll over repeats
if nargout >= 1
    nll_sd = sqrt(nansum(nll_var_vec));
end

if nargout > 2
    output.samples_used = samples_used;
    output.reps_used = reps_used;    
    output.funcalls = funcalls;
    output.allocate_reps = allocate_reps;
end

end