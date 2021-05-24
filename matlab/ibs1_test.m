%% test ibs2 one dataset
settings = get_model_settings('psycho') %'psycho', 'vstm','fourinarow'
dataset = load('data_psycho_s1')
stim = dataset.stim_all
% stim1 = cell2mat(stim)
stim1 = stim{1} %uncomment
resp = dataset.resp_all
% resp1 = cell2mat(resp)
resp1 = resp{1}; %uncomment
settings.Nsamples = 10; %this is for fixed sampling
[p_vec,Nreps,nll_exact,theta_inf,output_vec]=infer_theta('psycho','ibs_alloc_10',stim1,resp1,settings);
%% test all dataset
model = 'psycho';
%% test x point on ibs
method = "ibs_10";
method_split = split(method,"_");
xx = theta_inf
Ntrials= size(stim1,1);
thresh = settings.thresh*Ntrials;
if isempty(str2num(method_split{2})) %if it's ibs_alloc method, this condition should be true (empty)
    submethod = [method_split{1},'_', method_split{2}]
    Nsamples = str2num(method_split{3});
else 
    submethod = method_split{1}; %else, submethod = 'ibs'
    Nsamples = str2num(method_split{2});
end 
Nreps = Nsamples 
nll_ibs_arr = [];
nllvar_ibs_arr = [];
for i = 1:2000
    [nll_ibs,nll_sd_ibs,output_ibs]=estimate_nll_ibs(model,stim1,resp1,xx,Nreps,thresh)
    nll_ibs_arr = cat(1,nll_ibs_arr, nll_ibs);
    nllvar_ibs_arr = cat(1,nllvar_ibs_arr, nll_sd_ibs);

end 
nll_ibs = mean(nll_ibs_arr);
nllvar_ibs = mean(nllvar_ibs_arr);
%% test x point on ibs_alloc
method = "ibs_alloc_10";
method_split = split(method,"_");
if isempty(str2num(method_split{2})) %if it's ibs_alloc method, this condition should be true (empty)
    submethod = [method_split{1},'_', method_split{2}]
    Nsamples = str2num(method_split{3});
else 
    submethod = method_split{1}; %else, submethod = 'ibs'
    Nsamples = str2num(method_split{2});
end 
[nll_fix,~,output_fix] = estimate_nll_fixed1(model,stim1,resp1,xx,Nsamples); % use x0 to get vector p by calling fixed sampling
p_vec = output_fix.p_vec;

S_budget = round(sum(1./p_vec * Nsamples));%estimate the total number of sample (budget)
Nreps_alloc = round(S_budget * (1./sum(sqrt(dilog(p_vec)./p_vec)))*sqrt(p_vec.*dilog(p_vec))); %dilog(x) = Li_2(1-x);define nreps for each trial for ibs

nll_ibs2_arr = [];
nllvar_ibs2_arr = [];
for i = 1:2000
    [nll_ibs2,nll_sd_ibs2,output_ibs2]=estimate_nll_ibs(model,stim1,resp1,xx,Nreps_alloc,thresh)
    nll_ibs2_arr = cat(1,nll_ibs2_arr, nll_ibs2);
    nllvar_ibs2_arr = cat(1,nllvar_ibs2_arr, nll_sd_ibs2);

end 
nll_ibs2 = mean(nll_ibs2_arr);
nllvar_ibs2 = mean(nllvar_ibs2_arr);

%%
[nll_exact,~,~]=estimate_nll_exact(model,stim1,resp1,xx,1)
