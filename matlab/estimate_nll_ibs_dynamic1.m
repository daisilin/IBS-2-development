function [nll,nll_sd,output]=estimate_nll_ibs_dynamic1(model,stim,resp_real,theta,Nreps,thresh,p_initial, Nsamples,lookup_logp,dilog_p)
%ESTIMATE_NLL_IBS Negative log likelihood estimation via inverse binomial sampling.
% %thresh is for early stopping, if nll of data for each ibs run/repeat is smaller than a number; 
persistent samples_used;
persistent reps_used;
persistent funcalls;
persistent p_current; % current p vector from last function call, added for dynamic submethod
%persistent allocate_reps; % num of reps that should be allocated for this x
Ntrials = size(resp_real,1);
nll_trials = zeros(Ntrials,1);
K_tol = zeros(Ntrials,1);
if isempty(samples_used)
    samples_used = 0;
    reps_used = 0;
    funcalls = 0;
    p_current = p_initial;
    %allocate_reps = 0;
end
if nargin < 4 || isempty(theta)
    nll = []; nll_sd = []; 
    output.samples_used = samples_used;     
    output.reps_used = reps_used;
    output.funcalls = funcalls;
    output.p_current = p_current;
    samples_used = 0;
    reps_used = 0;
    funcalls = 0;
   % p_current = p_initial;
    return;
end

if nargin < 6 || isempty(thresh); thresh = Inf; end
if samples_used ~= 0
    S_budget = round(sum(1./p_current * Nsamples));%estimate the total number of sample (budget)
    dilog_lookup = interp1(lookup_logp,dilog_p,log(p_current));
    Nreps = round(S_budget * (1./sum(sqrt(dilog_lookup./p_current)))*sqrt(p_current.* dilog_lookup ));
else
    p_current = p_initial;
end
%Nreps=Nreps_dynamic_estimate(p_current, Nsamples);
%Nreps = ones(Ntrials,1)*3;
Nreps_max = max(Nreps);
nll_var_vec = zeros(1,Nreps_max); 
%nll_vec = zeros(1,Nreps_max); %use the largest repeat num; Nreps is a vector of length = N repeats,
if nargout > 1; nll_var_vec = zeros(1,Nreps_max); end

if isscalar(Nreps); Nreps = Nreps*ones(Ntrials,1); end

for iRep = 1:Nreps_max %loop over repeats
    resp = NaN(Ntrials,1);
    tries = zeros(Ntrials,1); % vector that counts how many tries/misses there have been before a match/hit for all trials; it gets reset for every repeat
    ind1 = iRep<=Nreps; %trials that need more reps = 1
    n_active = sum(ind1); %number of "active" trials (haven't matched) left (trials that needs repeat)
    ind_active=ind1;
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
         %nll_var_vec = compute_variance(iRep, K_active,ind_active,Nreps);
    end    
    samples_used = samples_used + (sum(tries)+n_active)/Ntrials;
    %nll_mat = cat(2,nll_mat,sum(nll_trials./K));
end

    %compute p_vec, added for dynamic submethod
    c = 0.0001;
    p_vec = (Nreps-0.5)./K_tol;% num hit (1) / num total tries(samples), assuming 0.5 hit 0.5 miss observed, a weaker prior
    p_current = c*p_current + (1-c)*p_vec;
funcalls = funcalls + 1;
reps_used = reps_used + mean(Nreps);
nll_mat = nll_trials./Nreps; % average nll over repeats
nll = sum(nll_mat); % sum over trials 
if nargout >= 1
    nll_sd = sqrt(nansum(nll_var_vec));
end

if nargout > 2
    output.samples_used = samples_used;
    output.reps_used = reps_used;    
    output.funcalls = funcalls;
    output.p_current = p_current;

    %output.allocate_reps = allocate_reps;

end

end