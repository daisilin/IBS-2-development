%% test ibs2 one dataset
settings = get_model_settings('psycho') %'psycho', 'vstm','fourinarow'
%dataset = load('data_psycho_s3')
dataset = load('data_psycho_s79') %model type

stim = dataset.stim_all
% stim1 = cell2mat(stim) 
stim1 = stim{9} %uncomment
resp = dataset.resp_all

% resp1 = cell2mat(resp)
resp1 = resp{9}; %uncomment

settings.Nsamples = 3; %this is for ibs
[p_vec,p0,tot_samples, nll_exact,theta_inf,output_vec]=infer_theta('psycho','ibs',stim1,resp1,settings,0.9);
%%
%% test ibs2 one simple bernoulli
b = 40;

settings.lb = [log(0.1) -2 0.01];
settings.ub = [log(10) 2 1];
settings.plb = [log(0.1) -1 0.01];
settings.pub = [log(5) 1 0.2];
settings.theta_real = [ ...
    [linspace(settings.plb(1),settings.pub(1),b)', 0.1*ones(b,1), 0.1*ones(b,1)]; ...
    [log(2)*ones(b,1), linspace(settings.plb(2),settings.pub(2),b)', 0.1*ones(b,1)]; ...
    [log(2)*ones(b,1), 0.1*ones(b,1), linspace(settings.plb(3),settings.pub(3),b)'] ...
    ];
settings.thresh = log(2);
settings.params = {'sigma','bias','lapse'};
settings.logflag = [1,0,0];
settings.methods = {'ibs','fixed','fixedb','exact'};
settings.samples{1} = [1 2 3 5 10 15 20 35 50];
settings.samples{2} = [1 2 3 5 10 15 20 35 50 100];
settings.samples{3} = [1 2 3 5 10 15 20 35 50 100];
settings.samples{4} = 0;
field1 = 'theta_real';  value1 = [0.9;0.9;0.1];
field2 = 'stim_all';  value2 = [1;1;1];
field3 = 'resp_all';  value3 = [1;1;1];

dataset = struct(field1,value1,field2,value2,field3,value3)

stim = dataset.stim_all
% stim1 = cell2mat(stim) 
% stim1 = stim{9} %uncomment
resp = dataset.resp_all
theta = dataset.theta_real
% resp1 = cell2mat(resp)
% resp1 = resp{9}; %uncomment
Ntrials = size(stim,1);
thresh = settings.thresh*Ntrials;
settings.Nsamples = 10; %this is for ibs
Nsamples=3;
p_true = [0.9;0.9;0.1];
p_wrong = [0.1;0.1;0.5];

lookup_logp = linspace(log(1e-7), log(1), 1e4); %sampled p points
dilog_p = dilog(exp(lookup_logp)); %value of dilog output
alpha=1;
highrep=20;
[nll_true,nll_sd_true,output_true]=estimate_nll_ibs2('bernoulli',stim,resp,theta,3,thresh,p_true, Nsamples,lookup_logp,dilog_p,highrep,alpha)
nll_var_true = nll_sd_true^2;
samples_true = output_true.samples_used;
[nll,nll_sd,output]=estimate_nll_ibs2('bernoulli',stim,resp,theta,3,thresh,p_wrong, Nsamples,lookup_logp,dilog_p,highrep,alpha)
nll_var = nll_sd^2;
samples = output.samples_used;
vs_ratio_empirical = (nll_var_true*samples_true)/(nll_var*samples);
var_ratio = nll_var_true/nll_var;
% [p_vec,p0,tot_samples, nll_exact,theta_inf,output_vec]=infer_theta('bernoulli','ibs_dynamic_10',stim,resp,settings,0.9);
%%
[nll,nll_sd,output]=estimate_nll_ibs_old('bernoulli',stim,resp,p_true,50,thresh)

%% using equation to get vs ratio
p0 = [0.9;0.9;0.1];
p_vec = [0.1;0.1;0.5];
vs_ratio= compute_vs_ratio(p0,p_vec);
%% get expected variance and expected samples used 
[exp_var,exp_sample] = get_var(p_vec,p0,Nsamples);
[exp_var0,exp_sample0] = get_var(p0,p0,Nsamples);
var_ratio_exp = exp_var0/exp_var;
%%
stim = dataset.stim_all;
resp = dataset.resp_all;

settings.Nsamples = 3; %this is for ibs
nll_exact_dy_arr = [];
nll_exact_st_arr = [];
nll_exact_ibs_arr = [];

nll_true_arr = [];

var_dy_arr = [];
var_st_arr = [];
var_ibs_arr = [];
for i = 1:10
    current_stim = stim{i};
    current_resp = resp{i};
    
    [p_vec_dy,tot_samples_dy, nll_exact_dy,theta_inf_dy,output_vec_dy]=infer_theta('vstm','ibs_dynamic_3',current_stim,current_resp,settings);
    nll_exact_dy_arr = cat(1,nll_exact_dy_arr,nll_exact_dy);
    var_dy_arr = cat(1,var_dy_arr,output_vec_dy(2));
    
    [p_vec_st,tot_samples_st, nll_exact_st,theta_inf_st,output_vec_st]=infer_theta('vstm','ibs_static_3',current_stim,current_resp,settings);
    nll_exact_st_arr = cat(1,nll_exact_st_arr,nll_exact_st);
    var_st_arr = cat(1,var_st_arr,output_vec_st(2));

    [p_vec,tot_samples, nll_exact_ibs,theta_inf,output_vec]=infer_theta('vstm','ibs_static_3',current_stim,current_resp,settings);
    nll_exact_ibs_arr = cat(1,nll_exact_ibs_arr,nll_exact_ibs);
    var_ibs_arr = cat(1,var_ibs_arr,output_vec(2));

    [~,~, nll_true,~,output_vec_exact]=infer_theta('vstm','exact',current_stim,current_resp,settings);
    nll_true_arr = cat(1,nll_true_arr,nll_true);

end 
%%
%%
stim = dataset.stim_all;
resp = dataset.resp_all;

settings.Nsamples = 3; %this is for ibs
nll_exact_dy_arr = [];
nll_exact_st_arr = [];
nll_exact_ibs_arr = [];

nll_true_arr = [];

var_dy_arr = [];
var_st_arr = [];
var_ibs_arr = [];
for i = 1:10
    current_stim = stim{i};
    current_resp = resp{i};
    
    [p_vec_dy,tot_samples_dy, nll_exact_dy,theta_inf_dy,output_vec_dy]=infer_theta('vstm','ibs_dynamic_3',current_stim,current_resp,settings);
    nll_exact_dy_arr = cat(1,nll_exact_dy_arr,nll_exact_dy);
    var_dy_arr = cat(1,var_dy_arr,output_vec_dy(2));
    
%     [p_vec_st,tot_samples_st, nll_exact_st,theta_inf_st,output_vec_st]=infer_theta('vstm','ibs_static_3',current_stim,current_resp,settings);
%     nll_exact_st_arr = cat(1,nll_exact_st_arr,nll_exact_st);
%     var_st_arr = cat(1,var_st_arr,output_vec_st(2));

%     [p_vec,tot_samples, nll_exact_ibs,theta_inf,output_vec]=infer_theta('vstm','ibs_3',current_stim,current_resp,settings);
%     nll_exact_ibs_arr = cat(1,nll_exact_ibs_arr,nll_exact_ibs);
%     var_ibs_arr = cat(1,var_ibs_arr,output_vec(2));

%     [~,~, nll_true,~,output_vec_exact]=infer_theta('vstm','exact',current_stim,current_resp,settings);
%     nll_true_arr = cat(1,nll_true_arr,nll_true);

end 
%% test all dataset
model = 'vstm';
%% test x point on ibs
method = "ibs_4";
method_split = split(method,"_");
xx = theta_inf; % the point to be evalutated 
%xx = [-2.2949,0.0850,0.0902];
%xx = [-1.00,0.7,0.9];
%xx = [-1.5,0.1,0.3];
%xx = [-0.8,0.2];
Ntrials= size(stim1,1);
thresh = settings.thresh*Ntrials;
if isempty(str2num(method_split{2})) %if it's ibs_alloc method, this condition should be true (empty)
    submethod = [method_split{1},'_', method_split{2}]
    Nsamples = str2num(method_split{3});
else 
    submethod = method_split{1}; %else, submethod = 'ibs'
    Nsamples = str2num(method_split{2});
end 
Nreps = Nsamples ;
nll_ibs_arr = [];
nllvar_ibs_arr = [];
samples_ibs_arr = [];
length = 500
for i = 1:length
    [nll_ibs,nll_sd_ibs,output_ibs]=estimate_nll_ibs(model,stim1,resp1,xx,Nreps,thresh);
    nll_ibs_arr = cat(1,nll_ibs_arr, nll_ibs);
    nllvar_ibs_arr = cat(1,nllvar_ibs_arr, nll_sd_ibs);
   % samples_ibs_arr =  cat(1,samples_ibs_arr, output_ibs.samples_used);
end 
nll_ibs = mean(nll_ibs_arr);
nllvar_ibs = sum(nllvar_ibs_arr);
samples_ibs =output_ibs.samples_used;
sem_ibs = std(nll_ibs_arr)/sqrt(length);
z_var_ibs = (nll_ibs - nll_ibs_arr)/std(nll_ibs_arr);
z_var_ibs_c = (nll_ibs - nll_ibs_arr)./nllvar_ibs_arr; % denominator-corrected zscore
var_per_sample_ibs= nllvar_ibs/samples_ibs;

%% test x point on ibs_alloc
method = "ibs_alloc_4";
method_split = split(method,"_");
if isempty(str2num(method_split{2})) %if it's ibs_alloc method, this condition should be true (empty)
    submethod = [method_split{1},'_', method_split{2}]
    Nsamples = str2num(method_split{3});
else 
    submethod = method_split{1}; %else, submethod = 'ibs'
    Nsamples = str2num(method_split{2});
end 
[nll_fix,~,output_fix] = estimate_nll_fixed1(model,stim1,resp1,xx,50); % use x0 to get vector p by calling fixed sampling
p_vec = output_fix.p_vec;

S_budget = round(sum(1./p_vec * Nsamples));%estimate the total number of sample (budget)
Nreps_alloc = round(S_budget * (1./sum(sqrt(dilog(p_vec)./p_vec)))*sqrt(p_vec.*dilog(p_vec))); %dilog(x) = Li_2(1-x);define nreps for each trial for ibs

nll_ibs2_arr = [];
nllvar_ibs2_arr = [];
for i = 1:length
    [nll_ibs2,nll_sd_ibs2,output_ibs2]=estimate_nll_ibs(model,stim1,resp1,xx,Nreps_alloc,thresh);
    nll_ibs2_arr = cat(1,nll_ibs2_arr, nll_ibs2);
    nllvar_ibs2_arr = cat(1,nllvar_ibs2_arr, nll_sd_ibs2);

end 
nll_ibs2 = mean(nll_ibs2_arr);
nllvar_ibs2 = sum(nllvar_ibs2_arr);
sem_ibs2 = std(nll_ibs2_arr)/sqrt(length);
samples_ibs2 =output_ibs2.samples_used;

z_var_ibs2 = (nll_ibs2 - nll_ibs2_arr)/std(nll_ibs2_arr);
z_var_ibs2_c = (nll_ibs2 - nll_ibs2_arr)./nllvar_ibs2_arr;

var_per_sample_ibs2= nllvar_ibs2/samples_ibs2;
%%
%% test x point on ibs_dynamic
method = "ibs_dynamic_4";
method_split = split(method,"_");
if isempty(str2num(method_split{2})) %if it's ibs_alloc method, this condition should be true (empty)
    submethod = [method_split{1},'_', method_split{2}]
    Nsamples = str2num(method_split{3});
else 
    submethod = method_split{1}; %else, submethod = 'ibs'
    Nsamples = str2num(method_split{2});
end 
[nll_fix,~,output_fix] = estimate_nll_fixed1(model,stim1,resp1,xx,50); % use x0 to get vector p by calling fixed sampling
p_initial = output_fix.p_vec;

S_budget = round(sum(1./p_initial * Nsamples));%estimate the total number of sample (budget)
Nreps_dynamic = round(S_budget * (1./sum(sqrt(dilog(p_initial)./p_initial)))*sqrt(p_initial.*dilog(p_initial))); %dilog(x) = Li_2(1-x);define nreps for each trial for ibs

nll_ibs3_arr = [];
nllvar_ibs3_arr = [];
%build loopup table for dilog calculation
lookup_logp = linspace(log(1e-6), log(1), 1e4); %sampled p points
dilog_p = dilog(exp(lookup_logp)); %value of dilog output
for i = 1:length
    [nll_ibs3,nll_sd_ibs3,output_ibs3]=estimate_nll_ibs_dynamic(model,stim1,resp1,xx,Nreps_dynamic,thresh,p_initial, Nsamples,lookup_logp,dilog_p);
    nll_ibs3_arr = cat(1,nll_ibs3_arr, nll_ibs3);
    nllvar_ibs3_arr = cat(1,nllvar_ibs3_arr, nll_sd_ibs3);

end 
nll_ibs3 = mean(nll_ibs3_arr);
nllvar_ibs3 = sum(nllvar_ibs3_arr);
sem_ibs3 = std(nll_ibs3_arr)/sqrt(length);
samples_ibs3 =output_ibs3.samples_used;

z_var_ibs3 = (nll_ibs3 - nll_ibs3_arr)/std(nll_ibs3_arr);
z_var_ibs3_c = (nll_ibs3 - nll_ibs3_arr)./nllvar_ibs3_arr;
var_per_sample_ibs3= nllvar_ibs3/samples_ibs3;

%%
[nll_exact,~,~]=estimate_nll_exact(model,stim1,resp1,xx,1)

%% plot

figure
qqplot(z_var_ibs_c);
title('QQ plot of zscore var ibs vs. standard normal');
figure
qqplot(z_var_ibs2_c);
title('QQ plot of zscore var ibs alloc vs. standard normal');
figure
qqplot(z_var_ibs3_c);
title('QQ plot of zscore var ibs dynamic vs. standard normal');
figure
histogram(z_var_ibs3_c);
