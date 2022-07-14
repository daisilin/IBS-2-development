function [p_vec,tot_samples,nll_exact,theta_inf,output_vec]=infer_theta(model,method,stim,resp,settings,alpha)
%INFER_THETA Optimization run for a single problem set.
% persistent x_00;
% if isempty(x_00)
%     x_00 = 0;
% end 
% Set default options for BADS optimization algorithm
method_split = split(method,"_");

badsopts = bads('defaults');
badsopts.UncertaintyHandling = ~strcmp(method,'exact'); %%changed method
badsopts.NoiseFinalSamples = 0;
badsopts.MaxFunEvals = 500;
if badsopts.UncertaintyHandling
    badsopts.SpecifyTargetNoise = true; %tell bads the noise of the point, which will improve the bad optimization 
end

% badsopts.NonlinearScaling = 'off';

% Extract settings
Nsamples = settings.Nsamples; %this is for both 'fixed' and 'ibs' method
% fprintf('nsample',Nsamples)
mult_hiprec = 10;
one_rep =1;
% switch method
%     case 'ibs'
% %         Nreps = Nsamples;
%         Ntrials = size(stim,1);
%         thresh = settings.thresh*Ntrials;
%         fun=@(x) estimate_nll_ibs(model,stim,resp,x,Nreps,thresh);
%         fun_hiprec=@(x) estimate_nll_ibs(model,stim,resp,x,Nreps*mult_hiprec,thresh);
%     case 'fixed'
%         fun = @(x) estimate_nll_fixed1(model,stim,resp,x,Nsamples);
%         fun_hiprec = @(x) estimate_nll_fixed1(model,stim,resp,x,Nsamples,mult_hiprec);
%     case 'fixedb'
%         fun = @(x) estimate_nll_fixed2(model,stim,resp,x,Nsamples);
%         fun_hiprec = @(x) estimate_nll_fixed2(model,stim,resp,x,Nsamples,mult_hiprec);
%     case 'exact'
%         fun=@(x) estimate_nll_exact(model,stim,resp,x,0);     
%         fun_hiprec=@(x) estimate_nll_exact(model,stim,resp,x,1);
%     otherwise
%         error(['Unknown method ''' method '''.']);
% end

% % Initialize estimation function if needed
% fun([]);

% Define grid of starting points
lb = settings.lb;
ub = settings.ub;
plb = settings.plb;
pub = settings.pub;
nvars = numel(pub);

switch nvars %each case correspond to different dimensions of parameters
    case 1; x0 = linspace(pb,pub,4)';
    case 2
        for j = 1:4
            x0(j,:)=plb+(pub-plb).*[1+mod(j-1,2),1+mod(floor((j-1)/2),2)]/3;
        end
    case 3
        for j = 1:8
            x0(j,:) = plb+(pub-plb).*[1+mod(j-1,2),1+mod(floor((j-1)/2),2),1+mod(floor((j-1)/4),2)]/3;
        end
end
% fprintf('x0',x0)
% method = {'ibs_alloc', 10};
if size(method_split,1) ==1
    submethod = method_split{1}        
elseif isempty(str2num(method_split{2})) %if it's ibs_static or dynamic method, this condition should be true (empty)
    submethod = [method_split{1},'_', method_split{2}]
    Nsamples = str2num(method_split{3});
else 
    submethod = method_split{1} %else, submethod = 'ibs'
    Nsamples = str2num(method_split{2});
end 
%method_split = strsplit(method);
%if length(method_split) > 1
%    Nsamples = str2num(method_split{2});
%end 


% Nreps = round(100 * sqrt(p_vec.*dilog(1-p_vec))); %define nreps for each trial for ibs

switch submethod
    case 'ibs'
        Nreps = Nsamples;
        Ntrials = size(stim,1);
        p_vec = ones(Ntrials,1);

        thresh = settings.thresh*Ntrials;
        fun=@(x) estimate_nll_ibs(model,stim,resp,x,Nreps,thresh);
        fun_hiprec=@(x) estimate_nll_ibs(model,stim,resp,x,Nreps*mult_hiprec,thresh);
    case 'ibs_static'
%       Nreps = Nsamples;
        Nreps_fixed = 100;
        Ntrials = size(stim,1);
        thresh = settings.thresh*Ntrials;
        [~,~,output_fix] = estimate_nll_fixed1(model,stim,resp,x0,Nreps_fixed); % use x0 to get vector p by calling fixed sampling
        p_vec = output_fix.p_vec;
        S_budget = round(sum(1./p_vec * Nsamples));%estimate the total number of sample (budget)
        Nreps = round(S_budget * (1./sum(sqrt(dilog(p_vec)./p_vec)))*sqrt(p_vec.*dilog(p_vec))); %dilog(x) = Li_2(1-x);define nreps for each trial for ibs
        Nreps(Nreps==0) = 1;
        fun=@(x) estimate_nll_ibs(model,stim,resp,x,Nreps,thresh);
        fun_hiprec=@(x) estimate_nll_ibs(model,stim,resp,x,Nreps*mult_hiprec,thresh);
    case 'ibs_dynamic'
        Nreps_fixed = 100;
        Ntrials = size(stim,1);
        thresh = settings.thresh*Ntrials;
        [~,~,output_fix] = estimate_nll_fixed1(model,stim,resp,x0,Nreps_fixed); % use x0 to get vector p by calling fixed sampling
        %build loopup table for dilog calculation
        lookup_logp = linspace(log(1e-7), log(1), 1e4); %sampled p points
        dilog_p = dilog(exp(lookup_logp)); %value of dilog output
        p_initial = output_fix.p_vec;
        S_budget = round(sum(1./p_initial * Nsamples));%estimate the total number of sample (budget)
        Nreps = round(S_budget * (1./sum(sqrt(dilog(p_initial)./p_initial)))*sqrt(p_initial.*dilog(p_initial))); %dilog(x) = Li_2(1-x);define nreps for each trial for ibs
        Nreps(Nreps==0) = 1;
        fun=@(x) estimate_nll_ibs2(model,stim,resp,x,Nreps,thresh,p_initial,Nsamples,lookup_logp,dilog_p,one_rep,alpha);
        fun_hiprec=@(x) estimate_nll_ibs2(model,stim,resp,x,Nreps*mult_hiprec,thresh,p_initial,Nsamples,lookup_logp,dilog_p,mult_hiprec,alpha);
    case 'fixed'
        fun = @(x) estimate_nll_fixed1(model,stim,resp,x,Nsamples);
        fun_hiprec = @(x) estimate_nll_fixed1(model,stim,resp,x,Nsamples,mult_hiprec);
    case 'fixedb'
        fun = @(x) estimate_nll_fixed2(model,stim,resp,x,Nsamples);
        fun_hiprec = @(x) estimate_nll_fixed2(model,stim,resp,x,Nsamples,mult_hiprec);
    case 'exact'
        fun=@(x) estimate_nll_exact(model,stim,resp,x,0);     
        fun_hiprec=@(x) estimate_nll_exact(model,stim,resp,x,1);
    otherwise
        error(['Unknown method ''' submethod '''.']);
end
% Initialize estimation function if needed
fun([]);

% Perform optimization from each starting point
Nopts = size(x0,1);
x_best = zeros(Nopts,nvars);
nLL = zeros(Nopts,1);
nLL_sd = zeros(Nopts,1);

% Nbench = 3;
% if Nbench > 0; bench_timing = bench(Nbench); end   % Measure computer speed
% t0 = tic;
for iOpt = 1:Nopts
    x_best(iOpt,:) = bads(fun,x0(iOpt,:),lb,ub,plb,pub,[],badsopts); %optimize with bads; fun is estimate_nll_ibs (defined above)
    % Evaluate candidate solution with higher precision
    [nLL(iOpt),nLL_sd(iOpt)] = fun_hiprec(x_best(iOpt,:));
end
% t_tot = toc(t0);
%if Nbench > 0; bench_timing = [bench_timing; bench(Nbench)]; end  % Measure computer speed

% Compute summed benchmark time (only numerical benchmarks)
% if Nbench > 0
%     bench_t = sum(sum(bench_timing(:,1:4)));
% else
%     bench_t = 0;
% end

% Get best solution
[~,idx_best] = min(nLL);
theta_inf = x_best(idx_best,:);
nLL_best = nLL(idx_best);
nLL_sd_best = nLL_sd(idx_best);

% Get function output (this also resets the function counters)
[~,~,output]=fun([]);

% Correct number of effective function calls to account for high prec
if strcmp(submethod,'ibs') || strcmp(submethod,'ibs_static') || strcmp(submethod,'ibs_dynamic') || strcmp(submethod,'fixed') || strcmp(submethod,'fixedb')
    output.funcalls = output.funcalls + Nopts*(mult_hiprec-1);
end
output_vec = [nLL_best,nLL_sd_best, output.samples_used,output.reps_used,output.funcalls];

%output_vec = [nLL_best,nLL_sd_best, output.samples_used,output.reps_used,output.funcalls,t_tot,bench_t];
tot_samples = output.samples_used;
if strcmp(submethod,'ibs_dynamic') 
    p_vec = output.p_current 
elseif strcmp(submethod,'ibs_static') 
    p_vec = p_vec;
% elseif strcmp(submethod,'ibs') 
%     p_vec = output.p_vec;
elseif strcmp(submethod,'exact') 
%     [nll,nll_sd,output]=estimate_nll_ibs(model,stim,resp,theta_inf,10);
    [~,~,output] = estimate_nll_fixed1(model,stim,resp,theta_inf,200); % use x0 to get vector p by calling fixed sampling
    p_vec = output.p_vec;
else
    p_vec = 0;
end 

[nll_exact,~,~]=estimate_nll_exact(model,stim,resp,theta_inf,1)
%nll_exact = nll_exact;
%nll_diff = nll_exact - output_vec(1);
%K_tot = output.K_tot;
%allocate_reps = output.allocate_reps;
end
