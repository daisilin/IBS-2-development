function [nll,nll_sd,output]=estimate_nll_ibs(model,stim,resp_real,theta,Nreps,thresh)
%ESTIMATE_NLL_IBS Negative log likelihood estimation via inverse binomial sampling.
% %thresh is for early stopping, if nll of data for each ibs run/repeat is smaller than a number; 
persistent samples_used;
persistent reps_used;
persistent funcalls;
%persistent K_tot; % trials x funcalls; nll for each trial, with all repeats averaged for each trial; each column is a different x (param setting)depending on bads
%persistent allocate_reps; % num of reps that should be allocated for this x

if isempty(samples_used)
    samples_used = 0;
    reps_used = 0;
    funcalls = 0;
    %K_tot = [];
    %allocate_reps = 0;
end
if nargin < 4 || isempty(theta)
    nll = []; nll_sd = []; 
    output.samples_used = samples_used;     
    output.reps_used = reps_used;
    output.funcalls = funcalls;
    %output.K_tot = K_tot;
    %output.allocate_reps = allocate_reps; 
    samples_used = 0;
    reps_used = 0;
    funcalls = 0;
    return;
end

Ntrials = size(resp_real,1);
nll_trials = zeros(Ntrials,1);
K_tol = zeros(Ntrials,1);

if nargin < 6 || isempty(thresh); thresh = Inf; end
if isscalar(Nreps); Nreps = Nreps*ones(Ntrials,1); end
Nreps_max = max(Nreps);
nll_var_vec = zeros(1,Nreps_max); 
%nll_vec = zeros(1,Nreps_max); %use the largest repeat num; Nreps is a vector of length = N repeats,
if nargout > 1; nll_var_vec = zeros(1,Nreps_max); end


for iRep = 1:Nreps_max %loop over repeats
    resp = NaN(Ntrials,1);
    tries = zeros(Ntrials,1); % vector that counts how many tries/misses there have been before a match/hit for all trials; it gets reset for every repeat
    ind1 = iRep<=Nreps; %trials that need more reps = 1
    ind_active=ind1;
    n_active = sum(ind1); %number of "active" trials (haven't matched) left (trials that needs repeat)
    n = sum(ind1);
    nll_sum = 0;
    nll_thresh = zeros(Ntrials,1);
    % Note: thresh here only makes sense in iterations that use all trials
    while n>0 && nll_sum<thresh % while there are trials that haven't been matched && nll of data for this repeat is smaller than a threshold
        resp(ind1) = feval(['generate_resp_' model],stim,theta,ind1); %generate synthetic responses from the model for active trials; 
        ind2 = ~(any(resp==resp_real,2)); %trial index where true and synthetic responses are not matched;  unmatched trials=1
        ind = any(ind1==1 & ind2==1,2); %index of unmatched and active trials (need repeat)
        ind1 = ind; % change ind1 to ind for resp(ind1) generating responses only for unmatched and active trials
        tries = tries+ind; %by adding ind, a logical vector where unmatched tirals = 1 ; we get a vector with Ntiral length, each idx counts how many misses for each (repeated) trial
        nll_trials(ind) = nll_trials(ind) + 1./(tries(ind)); %this accumulates over repeats; used for average calculation
        nll_thresh(ind) = nll_thresh(ind) + 1./(tries(ind)); %this gets zeroed for every repeat
        
        nll_sum = sum(nll_thresh);

        n = sum(ind); % count how many unmatched trials

    end
    K = tries+ind_active;
    K_tol = K_tol+K;
    K_active = K;
    K_active(K_active==0) = [];
    %Compute estimate of the variance if requested
    if nargout >= 1
        Ktab = -(psi(1,1:max(K_active(:)))' - psi(1,1));
        nll_var_vec(iRep) = nansum(Ktab(K_active)./Nreps(ind_active).^2);
    end      
    samples_used = samples_used + (sum(tries)+n_active)/Ntrials;
    %nll_mat = cat(2,nll_mat,sum(nll_trials./K));
end
    p_vec = (Nreps-0.5)./K_tol; % num hit (1) / num total tries(samples), assuming 0.5 hit 0.5 miss observed, a weaker prior

funcalls = funcalls + 1;
reps_used = reps_used + mean(Nreps);
nll_mat = nll_trials./Nreps;% average nll over repeats
nll = sum(nll_mat); 
if nargout >= 1
    nll_sd = sqrt(nansum(nll_var_vec));
end

if nargout > 2
    output.samples_used = samples_used;
    output.reps_used = reps_used;    
    output.funcalls = funcalls;
    output.p_vec = p_vec;
    %output.allocate_reps = allocate_reps;
end

end