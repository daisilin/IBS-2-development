function [nll,nll_sd,output]=estimate_nll_ibs(model,stim,resp_real,theta,Nreps,thresh)
%ESTIMATE_NLL_IBS Negative log likelihood estimation via inverse binomial sampling.
% %thresh is for early stopping, if nll of data for each ibs run/repeat is smaller than a number; 
persistent samples_used;
persistent reps_used;
persistent funcalls;
persistent K_tot; % trials x funcalls; nll for each trial, with all repeats averaged for each trial; each column is a different x (param setting)depending on bads
persistent allocate_reps; % num of reps that should be allocated for this x

if isempty(samples_used)
    samples_used = 0;
    reps_used = 0;
    funcalls = 0;
    K_tot = [];
    allocate_reps = 0;
end
if nargin < 4 || isempty(theta)
    nll = []; nll_sd = []; 
    output.samples_used = samples_used;
    output.reps_used = reps_used;
    output.funcalls = funcalls;
    output.K_tot = K_tot;
    output.allocate_reps = allocate_reps; 
    samples_used = 0;
    reps_used = 0;
    funcalls = 0;
    K_tot = [];
%     allocate_reps = 0; 
    return;
end

if nargin < 6 || isempty(thresh); thresh = Inf; end

% nll_vec = zeros(1,Nreps); % %a vector of nll for each ibs repeat
nll_vec = zeros(1,max(Nreps)); %use the largest repeat num; Nreps is a vector of length = Ntrial
if nargout > 1; nll_var_vec = zeros(1,max(Nreps)); end

nll_var_vec = zeros(1,max(Nreps)); 

Ntrials = size(resp_real,1);
K_reps = [];
K_tot = [];
% %initate variables for fix sampling
% fix_rep = 100; % 100 fix sampling 
% resp = NaN(Ntrials,1);
% fix_hits = zeros(Ntrials,1); %count misses for fix sampling 
% ind = true(Ntrials,1); 
% fix_ind = true(Ntrials,1); % logical where miss = 1 for fix sampling

% if size(K_tot,2)==1 % use fix sampling only when the first fun call (starting point of optimization)
%     for  fix_i = 1:fix_rep
        
%         resp(fix_ind) = feval(['generate_resp_' model],stim,theta,ind);%generate response for all trials rather than just miss trials
%         fix_ind(resp~=resp_real,2)=0; % logical; trial index with miss=0
%         fix_hits = fix_hits+fix_ind; % count how many hits for each trial
%     end 
%     p_0 = fix_hits/fix_rep; % a vector of likelihood, (hits/100), size = Ntrials x 1
%     allocate_reps = fix_rep * round(sqrt(p_0 .* dilog(1 - p_0))); %allocated repeats according to eq s14 in the paper 
% end 
% Loop over IBS repetitions (multiple IBS runs) for each trial
% Nreps = allocate_reps;
% for iTrial=1:Ntrials
for iRep=1:max(Nreps)
    resp = NaN(Ntrials,1);
    tries = zeros(Ntrials,1); % vector that counts how many tries/miss there have been before a match/hit for all trials; it gets reset for every repeat
%         tries_all = []; %count tries individual trial
    ind2= true(Ntrials,1); % initiate vector of trial idx where synthetic responses haven't matched resp_real yet
    ind1 = all(Nreps >= iRep); % logical vector of trial idx where sampling is required for this trial in this repeat
    ind = all(ind1==ind2,2); % both ind1 and ind2 are satisfied
    n = sum(ind); %number of trials (haven't matched) left

    while n>0 && nll_vec(iRep)<thresh % while there are trials that haven't been matched && nll of data for this repeat is smaller than a threshold
        resp(ind) = feval(['generate_resp_' model],stim,theta,ind); %generate synthetic response from the model; subsequent resps are for miss (and repeated) trials only
        ind = any(resp~=resp_real,2); %find trial index where true and synthetic responses are not matched; logical array where unmatched trials are set to 1
        tries = tries+ind; %by adding ind, a logical vector where unmatched tirals = 1 ; we get a vector with Ntiral length, each idx counts how many misses for each (repeated) trial
%             tries_all = tries_all + 1./(tries(ind)); % ntrials x 1; nll for each trial in this repeat by adding 1/misses
        nll_vec(iRep) = nll_vec(iRep) + sum(1./tries(ind)); % %get nll for this repeat, summed over unmatched trials; accumulate tries until all trial resp matched; 
        n=sum(ind);

    end
    %Compute estimate of the variance if requested
    if nargout > 1
        K = tries+1;
        Ktab = -(psi(1,1:max(K(:)))' - psi(1,1));
        nll_var_vec(iRep) = sum(Ktab(K));
    end

    %     for idx= 1:Ntrials % % calculate p vector!
    %         if tries(idx)+1 == 1
    %             k_val = 1;
    %         else 
    %             k_val = - 1./(tries(idx)+1);
    %         end 
    %         K_reps = cat(2,K_reps,k_val); % numtrials x numrepeats
    % 
    %     end 
    %     tries_inv = 1./(tries+1);
        samples_used = samples_used + mean(tries+1);
end 
% end 
% ave_K_reps = mean(K_reps,2); % average across repeats for each trial; numtrials x 1; (outside repeat for loop)
% if size(K_tot,2)==1
%     p_0 = exp(K_tot(:,1)); %convert loglike to likelihood p
%     allocate_reps = round(sqrt(p_0 .* dilog(1 - p_0))); %allocated repeats according to eq s14 in the paper 
% end 
funcalls = funcalls + 1;
reps_used = reps_used + max(Nreps);
% K_tot = cat(2,K_tot, ave_K_reps); %keep track of nll for all trials, and all x positions(all funcalls)

nll = mean(nll_vec); % return the average nll over repeats
if nargout >= 1
    nll_sd = sqrt(mean(nll_var_vec)/max(Nreps));
end

if nargout > 2
    output.samples_used = samples_used;
    output.reps_used = reps_used;    
    output.funcalls = funcalls;
    output.K_tot = K_tot;
    output.allocate_reps = allocate_reps;
end

end