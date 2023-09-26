function [nll,nll_var,output]=estimate_nll_aisp(func,stim,resp_real,theta,thresh,highrep,Nreps)
%ESTIMATE_NLL_ONLY only do estimation for Negative log likelihood
%estimation via inverse binomial sampling, without calculating allocations
%model is the generative model
%thresh is for early stopping
%Nreps can be a vector (in adaptive methods), or a scalar (in vanilla IBS)

%initialize variables
Ntrials = size(resp_real,1);
nll_trials = zeros(Ntrials,1);
K_tol = zeros(Ntrials,1);
samples_current = 0;
persistent samples_used;
if isempty(samples_used)
    samples_used = 0;
end 
    
%set thresh to infinity
if nargin < 6 || isempty(thresh); thresh = Inf; end
%for adaptive allocation, use the highest repeat number

%calculate Nreps
if isscalar(Nreps); Nreps = Nreps*ones(Ntrials,1);end
Nreps = Nreps*highrep;
Nreps_max = max(Nreps);
K_mtx = ones(Ntrials,Nreps_max);

%do IBS Nreps_max times
for iRep = 1:Nreps_max 
    resp = NaN(Ntrials,1);
    tries = zeros(Ntrials,1); % vector that counts how many misses there have been before a hit for all trials; it gets reset for every repeat
    ind1 = iRep <= Nreps; %index vector; True for active trials that need more repeats, False for trials that are done with thier max repeat
    ind_active = ind1; % variable for active trials, since ind1 will be changed in the loop
    n_active = sum(ind1); %number of "active" trials (no hit) left (trials that still need repeat)
    n = sum(ind1);
    nll_sum = 0;
    nll_thresh = zeros(Ntrials,1);
    
    % loop while there are trials that haven't been matched && nll of data for this repeat is smaller than a threshold
    % Note: thresh here only makes sense in iterations that use all trials
    while n>0 && nll_sum<thresh 
%         resp(ind1) = ibs_fun(data,theta,model,ind1); %simulate responses from the model for active trials && no-hit trials only; 
        resp(ind1) = func(ind1); %simulate responses from the model for active trials && no-hit trials only; 

        ind2 = any(resp==resp_real,2); % index vector for trials that has no hit;  no hit trials = 0
        ind = any(ind1==1 & ind2==0,2); % index vector for no-hit && active trials
        ind1 = ind; % ind1 is now active && no-hit trial indices, so that resp(ind1) generates responses for no-hit && active trials
        tries = tries+ind; % each idx adds the number of miss for active trials only; a vector of Ntrials length
        nll_trials(ind) = nll_trials(ind) + 1./(tries(ind)); %this nll vector accumulates over repeats; use for average calculation
        nll_thresh(ind) = nll_thresh(ind) + 1./(tries(ind)); %this nll vector gets zeroed for every repeat; use for nll_sum calculation
        nll_sum = sum(nll_thresh);
        n = sum(ind); % count how many no-hit && active trials

    end
    K = tries+ind_active; % K = misses + 1. +1 for active trials
    K_nan = K;
    K_nan(K_nan==0)=1; 
    K_mtx(:,iRep) = K_nan; 
    K_tol = K_tol+K;

    samples_used = samples_used + (sum(tries)+n_active)/Ntrials;
    samples_current = samples_current+(sum(tries)+n_active)/Ntrials;
end
p_vec = Nreps./K_tol;
    %Compute estimate of the variance if requested
    if nargout >= 1
         Ktab = -(psi(1,1:max(K_mtx(:)))' - psi(1,1));
         LLvar = Ktab(K_mtx);   
         nlogLvar = sum(LLvar,2)./Nreps.^2; %variance for each trial,summed over repeats
         nll_var = sum(nlogLvar);
         var_per_sample = nll_var/samples_current;
    end  
% average nll over repeats
nll_mat = nll_trials./Nreps;
nll = sum(nll_mat); 
if nargout > 2
    output.samples_used = samples_used;
    output.samples_current = samples_current;
    output.var_per_sample = var_per_sample;
    output.p_vec = p_vec;
    output.nlogLvar_trials = nlogLvar;
    output.nll_mat = nll_mat;
end