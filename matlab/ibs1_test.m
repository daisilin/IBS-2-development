%% test ibs2 one dataset
settings = get_model_settings('psycho') %'psycho', 'vstm','fourinarow'
dataset = load('data_psycho_s1')
stim = dataset.stim_all
% stim1 = cell2mat(stim)
stim1 = stim{1} %uncomment
resp = dataset.resp_all
% resp1 = cell2mat(resp)
resp1 = resp{1}; %uncomment
settings.Nsamples = 50; %this is for fixed sampling
[p_vec,Nreps,K_tot,allocate_reps,theta_inf,output_vec]=infer_theta('psycho',{'ibs_alloc',50},stim1,resp1,settings);
%% test all dataset
recover_theta('psycho','ibs2',1,5,100)