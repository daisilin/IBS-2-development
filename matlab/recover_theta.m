function recover_theta(model,method,alpha,proc_id,Nsamples,Ndatasets)
%RECOVER_THETA Run batch parameter recovery for given set parameters.

if nargin < 4; Nsamples = []; end
if nargin < 5; Ndatasets = []; end

method_split = split(method,"_");
settings = get_model_settings(model);
if method_split == "exact"
    Nsamples = 1; 
elseif isempty(str2num(method_split{2})) 
    %submethod = [method_split{1},'_', method_split{2}]
    Nsamples = str2num(method_split{3});
else  
    %submethod = method_split{1};
    Nsamples = str2num(method_split{2});
end 
settings.Nsamples = Nsamples;

% Add required folders

mypath = fileparts(mfilename('fullpath'));
addpath([mypath filesep 'bads']);
addpath([mypath filesep 'datasets']);
if strcmp(model,'vstm'); addpath([mypath filesep 'CircStat2012a']); end

% Loop over iterations (add case where method is not string...)
% methodtype = isa(method,'char') 

% if methodtype == 1
%     theta_filename = ['theta_' model '_' method '_' num2str(proc_id) '.txt'];
%     output_filename = ['output_' model '_' method '_' num2str(proc_id) '.txt'];
% else
%theta_filename = ['theta_' model '_' cell2mat(method_split) '_' num2str(proc_id) '.txt'];
theta_filename = ['theta_' model '_' method '_' num2str(proc_id) '.txt'];

%output_filename = ['output_' model '_' cell2mat(method_split) '_' num2str(proc_id) '.txt'];
output_filename = ['output_' model '_' method '_' num2str(proc_id) '.txt'];
%nll_exact_filename = ['nll_exact_' model '_' cell2mat(method_split) '_' num2str(proc_id) '.txt'];
nll_exact_filename = ['nll_exact_' model '_' method '_' num2str(proc_id) '.txt'];

%pvec_filename = ['p_vec_' model '_' cell2mat(method_split) '_' num2str(proc_id) '.txt'];
%Nreps_filename = ['Nreps' model '_' cell2mat(method_split) '_' num2str(proc_id) '.txt'];
pvec_filename = ['p_vec_' model '_' method '_' num2str(proc_id) '.txt'];
tot_samples_filename = ['tot_samples_' model '_' method '_' num2str(proc_id) '.txt'];
%nll_diff_filename = ['nll_diff' model '_' cell2mat(method_split) '_' num2str(proc_id) '.txt'];

%save('/scratch/xl1005/IBS-2-development/results/theta_filename.txt','theta_filename')
%save('/scratch/xl1005/IBS-2-development/results/output_filename.txt','output_filename')
%save('/scratch/xl1005/IBS-2-development/results/xbest_filename.txt','xbest_filename')
%save('/scratch/xl1005/IBS-2-development/results/pvec_filename.txt','pvec_filename')
%save('/scratch/xl1005/IBS-2-development/results/Nreps_filename.txt','Nreps_filename')

% If output file already exist, continue from previous run
theta_inf = []; output_vec = [];nll_exact = []; p_vec=[];tot_samples =[];
if exist(theta_filename,'file') && exist(output_filename,'file')
    try
        theta_inf = dlmread(theta_filename);
        output_vec = dlmread(output_filename);
        nll_exact = dlmread(nll_exact_filename);
        %nll_diff = dlmread(nll_diff_filename);
        p_vec = dlmread(pvec_filename);
        tot_samples = dlmread(tot_samples_filename);
    catch
        % File(s) corrupted, need to restart from scratch
    end
end

% Loop over simulated datasets
iStart = size(theta_inf,1)+1;

% Load fake datasets for given model and parameter setting
datafile = [mypath filesep 'datasets' filesep 'data_' model '_s' num2str(proc_id) '.mat'];
data = load(datafile);
if isempty(Ndatasets); Ndatasets = numel(data.stim_all); end

for i=iStart:Ndatasets
    rng(i);
    stim=data.stim_all{i};
    resp=data.resp_all{i};
    [p_vec(i,:),tot_samples(i,:),nll_exact(i,:),theta_inf(i,:),output_vec(i,:)] = ...
        infer_theta(model,method,stim,resp,settings,alpha);
    dlmwrite(theta_filename,theta_inf,'Delimiter','\t')
    dlmwrite(output_filename,output_vec,'Delimiter','\t')
    dlmwrite(nll_exact_filename,nll_exact,'Delimiter','\t')
   % dlmwrite(nll_diff_filename,nll_diff,'Delimiter','\t')

    dlmwrite(pvec_filename,p_vec,'Delimiter','\t')
    dlmwrite(tot_samples_filename,tot_samples,'Delimiter','\t')
end
  
end


