%% use exact params
folder = '/Users/daisy/Desktop/projects/ibs2/ibs-development-master/results/psycho/exact';
%folder = '/Users/daisy/Desktop/projects/ibs2/ibs-development-master/results/vstm/exact';
% Check to make sure that folder actually exists.  Warn user if it doesn't.
if ~isfolder(folder) &&  ~isfolder(folder)
  errorMessage = sprintf('Error: The following folder does not exist:\n%s', myFolder);
  uiwait(warndlg(errorMessage));
  return;
end
% Get a list of all files in the folder with the desired file name pattern.
file_pattern = fullfile(folder, 'theta*.txt'); % Change to whatever pattern you need.
files = dir(file_pattern);
thetas = [];
%loop through all (80 or 100 for vstm and psycho) datasets and get thetas 
for k = 1 : length(files)
    filename = files(k).name;
    fullname = fullfile(folder, filename);
    theta_current_all = importdata(fullname);
    step = 2;
    %loop through 100 (or selected subset) tasks in each dataset
    for j = 1:step:100
        theta_current = theta_current_all(j,:); %take the first param set
        thetas = cat(1,thetas,theta_current);
    end
end

%build lookup table for Nreps calculation later
lookup_logp = linspace(log(1e-7), log(1), 1e4); %sampled p points
dilog_p = dilog(exp(lookup_logp)); %value of dilog output
%% plots for x exact
all_nlls_best = [];
model = 'psycho';
settings = get_model_settings(model);
%loop through params for each dataset
task_num = 0;
for n = 1:length(files)
    
    dataset_name = ['data_' model '_' 's' num2str(n)];
    ibs_name = 'ibs_50';
    nsamples = 50;
    task_num = 100;     
    dataset = load(dataset_name)
    %loop through 100 tasks in each dataset
    
    for m = 1:step:task_num
        task_num = task_num + 1;
        stim_num = m;
        stim = dataset.stim_all;
        stim1 = stim{stim_num};
    
        resp = dataset.resp_all;
        resp1 = resp{stim_num}; 
        Ntrials = size(resp1,1);
        thresh = settings.thresh*Ntrials;
        settings.Nsamples = nsamples; %this is for ibs
        %get exact params
        theta = thetas(task_num,:);
        %theta = dataset.theta_real;
        [~,~,output_best]=estimate_nll_test(model,stim1,resp1,theta,thresh,1,5);
        nll_trials_best = output_best.nll_mat;
        ll_trials_best = -(nll_trials_best);
        all_nlls_best = cat(1,all_nlls_best,ll_trials_best);    
    end
    
end
% convert ll vector to likelihoods and use s16 to calculate the ratio, x best
%to store all variance ratios
ratios_xbest = [];
%set trial number
if strcmp('vstm',model) 
    trial_per_data = 400;
else
    trial_per_data = 600;
end
total_trial_num = length(all_nlls_best);
%get ratio for each dataset
for dataset_start =  1:trial_per_data:total_trial_num
    current_trials = all_nlls_best(dataset_start:dataset_start+trial_per_data-1);
    like_trials_xbest = exp(current_trials);
    v_ratio_xbest= compute_vratio_s16(like_trials_xbest);
    ratios_xbest = cat(1,ratios_xbest,v_ratio_xbest);
end

%plot histogram
histogram(ratios_xbest,10)
xlabel('variance ratio','FontSize', 10)
ylabel('frequency','FontSize', 8)
xlim([0.2 1])
% ylim([0 400])
figtitle = ['histogram of variance ratio,' model];
sgtitle(figtitle,'FontSize', 10)
%%
%% plots for variance ratio inital vs best
all_vars_best = [];
all_vars_init = [];
model = 'psycho';

settings = get_model_settings(model);
% get x0
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
%loop through params for each dataset
task_num = 0;
if strcmp('vstm',model)
    init_total = 4;
end
if strcmp('psycho',model)
    init_total = 8;
end

for init_id = 1:2:init_total
    theta_init = x0(init_id,:); 
    for n = 1:2:length(files)

        dataset_name = ['data_' model '_' 's' num2str(n)];
        ibs_name = 'ibs_5';
        nsamples = 5;
        task_num = 100;     
        dataset = load(dataset_name);
        %loop through 100 tasks in each dataset

        for m = 1:1:task_num
            task_num = task_num + 1;
            stim_num = m;
            stim = dataset.stim_all;
            stim1 = stim{stim_num};

            resp = dataset.resp_all;
            resp1 = resp{stim_num}; 
            Ntrials = size(resp1,1);
            thresh = settings.thresh*Ntrials;
            settings.Nsamples = nsamples; %this is for ibs
            %get exact params
            %theta = thetas(task_num,:);
            %theta_init = x0(1,:);
            theta = dataset.theta_real;
            %estimate for p vectors 
            %[~,~,output_best]=estimate_nll_test(model,stim1,resp1,theta,thresh,1,nsamples);
            if strcmp('vstm',model)
                [exact_ll,exact_ll_trials] = compute_nll_vstm(stim1,resp1,theta);
                pvec_best = exp(exact_ll_trials);
                [exact_init_ll,exact_init_ll_trials] = compute_nll_vstm(stim1,resp1,theta_init);
                pvec_init = exp(exact_init_ll_trials);
            end 
            if strcmp('psycho',model)
                [exact_ll,exact_ll_trials] = compute_nll_psycho(stim1,resp1,theta);
                pvec_best = exp(exact_ll_trials);
                [exact_init_ll,exact_init_ll_trials] = compute_nll_psycho(stim1,resp1,theta_init);
                pvec_init = exp(exact_init_ll_trials);
            end    
            %[~,~,output_init]=estimate_nll_test(model,stim1,resp1,theta_init,thresh,1,nsamples);
            %use fixed sampling to get init pvec
            %[~,~,output_init] = estimate_nll_fixed1(model,stim1,resp1,theta_init,50);
            % use p vectors to allocate resources
            %pvec_best = output_best.p_vec;
            [Nreps_best,budget_best] = estimate_Nreps(nsamples,lookup_logp,dilog_p,pvec_best,true,pvec_best);
%             pvec_init = output_init.p_vec;
            [Nreps_init,budget_init] = estimate_Nreps(nsamples,lookup_logp,dilog_p,pvec_init,true,pvec_best);

            [~,var_best,output_alloc_best]=estimate_nll_test(model,stim1,resp1,theta,thresh,1,Nreps_best);
            [~,var_init,output_alloc_init]=estimate_nll_test(model,stim1,resp1,theta,thresh,1,Nreps_init);
            nll_trials_best = output_alloc_best.nll_mat;
            nll_trials_init = output_alloc_init.nll_mat;

            all_vars_best = cat(1,all_vars_best,var_best);
            all_vars_init = cat(1,all_vars_init,var_init);
        end

    end
end
% convert ll vector to likelihoods and use s16 to calculate the ratio, x best
%to store all variance ratios
ratios_init_vs_best = [];
%set trial number
if strcmp('vstm',model) 
    trial_per_data = 400;
else
    trial_per_data = 600;
end
total_trial_num = length(all_vars_best);
%get ratio for each dataset
ratios_init_vs_best = all_vars_best./all_vars_init;
%plot histogram
histogram(ratios_init_vs_best,15)

edges = [0.6:0.02:1.15];
histogram(ratios_init_vs_best,edges);
xlabel('variance ratio','FontSize', 10)
ylabel('frequency','FontSize', 8)
%xlim([0 2])
% ylim([0 400])
figtitle = ['histogram of variance ratio,' model];
sgtitle(figtitle,'FontSize', 10)

%% AISP
% Get a list of all files in the folder with the desired file name pattern.
folder = '/Users/daisy/Desktop/projects/ibs2/AISP-master/categorizationWill/pars';
file_pattern = fullfile(folder, 'pars_ibs_Bayes_*.mat'); % Change to whatever pattern you need.
files = dir(file_pattern);
thetas = [];
liks = [];
% dataset_name = ['/Users/daisy/Desktop/projects/ibs2/AISP-master/categorizationWill/pars/pars_ibs_Bayes_' num2str(sub) '_' num2str(rep)];
% data = load(dataset_name);
%loop through all datasets and get thetas and nlls
for k = 1 : length(files)
    filename = files(k).name;
    fullname = fullfile(folder, filename);
    theta_current_all = importdata(fullname);

    theta_current = theta_current_all.pars; %take the param set
    thetas = cat(1,thetas,theta_current);
    lik_current = theta_current_all.likelihood; %take the param set
    liks = cat(1,liks,lik_current);
end
%get ratio for each dataset
ratios_xbest = [];
for dataset_start =  1:length(files)

    v_ratio_xbest= compute_vratio_s16(liks);
    ratios_xbest = cat(1,ratios_xbest,v_ratio_xbest);
end