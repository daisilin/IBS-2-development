function [nll,nll_var,output]=estimate_for_all(model,stim,resp_real,theta,Nsamples,thresh,lookup_logp,dilog_p,highrep,alpha,p_initial)
%estimation of negative log likelihood for all IBS methods including
%dynamic allocations
%Nsamples is number of IBS repeats. 
%p_initial is a vector of likelihoods for each trial (estimated through
%fixed sampling); p_initial is empty for IBS vanilla and static

% persistent samples_used;
persistent reps_used;
persistent funcalls;
persistent p_current;

if isempty(funcalls)
%     samples_used = 0;
    reps_used = 0;
    funcalls = 0;
end
if nargin < 4 || isempty(theta)
    nll = []; nll_var = []; 
    output.reps_used = reps_used;
    output.funcalls = funcalls;
    output.p_current = p_current;
    reps_used = 0;
    funcalls = 0;
    return;
end
% for vanilla IBS and static, p_current=1
if ~exist('p_initial','var') || isempty(p_initial)
    Nreps = Nsamples;
    alpha = 1;
    p_current = 1;
elseif funcalls ~= 0 %if there is p_initial input (it's dynamic method) and the function has been called
    Nreps = estimate_Nreps(Nsamples,lookup_logp,dilog_p,p_current)
else %if there is p_initial and it's the first function call, use the S_budget input directly
    p_current = p_initial;
    Nreps = estimate_Nreps(Nsamples,lookup_logp,dilog_p,p_current)
end

if model == "bernoulli"; theta = p_initial;end

if nargin < 6 || isempty(thresh); thresh = Inf; end

[nll,nll_var,output_only]= estimate_nll_only(model,stim,resp_real,theta,thresh,lookup_logp,dilog_p,highrep,p_current,Nreps);

% update p_current if it's dynamic method
if p_current ~=1; p_current = alpha*p_current + (1-alpha)*p_vec; end

funcalls = funcalls + 1;
reps_used = reps_used + mean(Nreps);

if nargout > 2
    output.samples_used = output_only.samples_used;
    output.reps_used = reps_used;    
    output.funcalls = funcalls;
    output.p_current = p_current;
%     output.nlogLvar_trials = nlogLvar;
end

end