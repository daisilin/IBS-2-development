function [nll,nll_var,output]=estimate_nll_ibs2(model,stim,resp_real,theta,reps,thresh,p_initial, Nsamples,lookup_logp,dilog_p,highrep,alpha)
%ESTIMATE_NLL_IBS Negative log likelihood estimation via inverse binomial sampling.
% %thresh is for early stopping, if nll of data for each ibs run/repeat is smaller than a number; 
persistent samples_used;
persistent reps_used;
persistent funcalls;
persistent p_current;

%persistent K_tot; % trials x funcalls; nll for each trial, with all repeats averaged for each trial; each column is a different x (param setting)depending on bads

if isempty(samples_used)
    samples_used = 0;
    reps_used = 0;
    funcalls = 0;
end
 
if nargin < 4 || isempty(theta)
    nll = []; nll_var = []; 
    output.samples_used = samples_used;     
    output.reps_used = reps_used;
    output.funcalls = funcalls;
    output.p_current = p_current;
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
total_samples = 0;
% 
if ~exist('p_initial','var') || isempty(p_initial)
    p_initial = ones(Ntrials,1);
    
end


if p_initial == 1
    Nreps = reps;
    alpha = 1;
    p_current = 1;
elseif samples_used ~= 0
    S_budget = round(sum(1./p_current * Nsamples));%estimate the total number of sample (budget)
    dilog_lookup = interp1(lookup_logp,dilog_p,log(p_current));
    Nreps = round(S_budget * (1./sum(sqrt(dilog_lookup./p_current)))*sqrt(p_current.* dilog_lookup ));
     %Nreps = round(S_budget * (1./sum(sqrt(dilog(p_current)./p_current)))*sqrt(p_current.*dilog(p_current))); %dilog(x) = Li_2(1-x);define nreps for each trial for ibs
    Nreps(Nreps==0) = 1;
else
    p_current = p_initial;
    S_budget = round(sum(1./p_current * Nsamples));%estimate the total number of sample (budget)
    dilog_lookup = interp1(lookup_logp,dilog_p,log(p_current));
    Nreps = round(S_budget * (1./sum(sqrt(dilog_lookup./p_current)))*sqrt(p_current.* dilog_lookup ));
    Nreps(Nreps==0) = 1;
%     Nreps=reps;
end
Nreps = Nreps*highrep;

if model == "bernoulli"; theta = p_initial;end

if nargin < 6 || isempty(thresh); thresh = Inf; end
% if isscalar(reps); Nreps = reps*ones(Ntrials,1); end % NEED TO FIX, PASS
% EVERYTHING AS VECTOR, NO SCALAR INTO THE FUNCTION

Nreps_max = max(Nreps);
% nll_var_vec = zeros(1,Nreps_max); 
K_mtx = ones(Ntrials,Nreps_max);

% if nargout > 1; nll_var_vec = zeros(1,Nreps_max); end


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
%         ind2 = ~(any(resp==resp_real,2)); %trial index where true and synthetic responses are not matched;  unmatched trials=1
        ind2 = any(resp==resp_real,2); %trial index where true and synthetic responses are not matched;  unmatched trials=0

        ind = any(ind1==1 & ind2==0,2); %index of unmatched and active trials (need repeat)
        ind1 = ind; % change ind1 to ind for resp(ind1) generating responses only for unmatched and active trials
        tries = tries+ind; %by adding ind, a logical vector where unmatched tirals = 1 ; we get a vector with Ntiral length, each idx counts how many misses for each (repeated) trial
        nll_trials(ind) = nll_trials(ind) + 1./(tries(ind)); %this accumulates over repeats; used for average calculation
        nll_thresh(ind) = nll_thresh(ind) + 1./(tries(ind)); %this gets zeroed for every repeat
        
        nll_sum = sum(nll_thresh);

        n = sum(ind); % count how many unmatched trials

    end
    K = tries+ind_active;
    K_nan = K;
    K_nan(K_nan==0)=1;
    K_mtx(:,iRep) = K_nan; 
    K_tol = K_tol+K;
    K_active = K;
    K_active(K_active==0) = [];

    total_samples = total_samples + (sum(tries)+n_active)/Ntrials;
    samples_used = samples_used + (sum(tries)+n_active)/Ntrials;
    %nll_mat = cat(2,nll_mat,sum(nll_trials./K));
end
    %Compute estimate of the variance if requested
    if nargout >= 1
         Ktab = -(psi(1,1:max(K_mtx(:)))' - psi(1,1));
         LLvar = Ktab(K_mtx);   
         nlogLvar = sum(LLvar,2)./Nreps.^2; %var per trial,summed over repeats
%          nll_var_vec_i = compute_variance(K_active,ind_active,Nreps);
%          nll_var_vec(iRep) = nll_var_vec_i;
    end  
%     alpha = 0.5;
%     p_vec = (Nreps-0.5)./K_tol;% num hit (1) / num total tries(samples), assuming 0.5 hit 0.5 miss observed, a weaker prior
    p_vec = Nreps./K_tol;
    if p_current ~=1; p_current = alpha*p_current + (1-alpha)*p_vec; end

%     p_current = alpha*p_current + (1-alpha)*p_vec;
funcalls = funcalls + 1;
reps_used = reps_used + mean(Nreps);
nll_mat = nll_trials./Nreps;% average nll over repeats
nll = sum(nll_mat); 
if nargout >= 1
    %nll_sd = nanmean(sqrt(nll_var_vec));
%     nll_sd = sqrt(nanmean(nll_var_vec));
      nll_var = sum(nlogLvar);
%     nll_sd = sqrt(nanmean(nll_var_vec)/mean(Nreps));

end

if nargout > 2
    output.samples_used = samples_used;
    output.reps_used = reps_used;    
    output.funcalls = funcalls;
    output.p_current = p_current;
    output.nlogLvar_trials = nlogLvar;
end

end