P = '/Users/daisy/Desktop/projects/ibs2/results2/psycho/50/';
F = sprintf('*%i*.txt','nll_exact_psycho_ibs_alloc');
%S = dir(fullfile(P,F));
%N = S.name;
%I = imread(fullfile(P,N));

S = dir(fullfile(P,'*.txt'));
%N = {S.name};

%X = cell(1,100);
I = cell(1,720);

for i = 1:720
    N = {S.name};
    X = ~cellfun('isempty',strfind(N,'nll_exact_psycho_ibs_alloc'));

    I{i} = dlmread(fullfile(P,N{X}));
end 
%%
ibs_alloc = dlmread('/Users/daisy/Desktop/projects/ibs2/results2/psycho/100/nll_exact_psycho_ibs_alloc50_50.txt')

ibs = dlmread('/Users/daisy/Desktop/projects/ibs2/results2/psycho/ibs100/nll_exact_psycho_ibs_50.txt')

diff = ibs_alloc - ibs(1:38)
%%
numfiles = 120;
vstm_ibs_5= cell(1, 80);
vstm_ibs_alloc_5= cell(1, 80);
psy_ibs_5= cell(1, numfiles);
psy_ibs_alloc_5= cell(1, numfiles);
psy_ibs_10= cell(1, numfiles);
psy_ibs_alloc_10= cell(1, numfiles);
%vstm_ibs_5= cell(1, numfiles);

for k = 1:numfiles
  %file_vstm_alloc_5 = sprintf('nll_exact_vstm_ibs_alloc_5_%d.txt', k);
  %file_four_alloc_5 = sprintf('nll_exact_fourinarow_ibs_alloc_5_%d.txt', k);
  file_psy_alloc_5 = sprintf('nll_exact_psycho_ibs_alloc_5_%d.txt', k);
  file_psy_ibs_5 = sprintf('nll_exact_psycho_ibs_5_%d.txt', k);
  file_psy_alloc_10 = sprintf('nll_exact_psycho_ibs_alloc_10_%d.txt', k);
  file_psy_ibs_10 = sprintf('nll_exact_psycho_ibs_10_%d.txt', k);
  %file_vstm_ibs_5 = sprintf('nll_exact_vstm_ibs_5_%d.txt', k);
 %file_four_ibs_5 = sprintf('nll_exact_fourinarow_ibs_5_%d.txt', k);
  %vstm_ibs_5{k} = importdata(file_vstm_ibs_5);
  %vstm_ibs_alloc_5{k} = importdata(file_vstm_alloc_5);
  %four_ibs_5{k} = importdata(file_four_ibs_5);
  %four_ibs_alloc_5{k} = importdata(file_four_alloc_5);
  psy_ibs_5{k} = importdata(file_psy_ibs_5);
  psy_ibs_alloc_5{k} = importdata(file_psy_alloc_5);
  psy_ibs_10{k} = importdata(file_psy_ibs_10);
  psy_ibs_alloc_10{k} = importdata(file_psy_alloc_10);
end

for j = 1:80
  file_vstm_alloc_5 = sprintf('nll_exact_vstm_ibs_alloc_5_%d.txt', j);
  %file_four_alloc_5 = sprintf('nll_exact_fourinarow_ibs_alloc_5_%d.txt', k);
  %file_psy_alloc_5 = sprintf('nll_exact_psycho_ibs_alloc_5_%d.txt', k);

  file_vstm_ibs_5 = sprintf('nll_exact_vstm_ibs_5_%d.txt', j);
 %file_four_ibs_5 = sprintf('nll_exact_fourinarow_ibs_5_%d.txt', k);
  %file_psy_ibs_5 = sprintf('nll_exact_psycho_ibs_5_%d.txt', k);

  vstm_ibs_5{j} = importdata(file_vstm_ibs_5);
  vstm_ibs_alloc_5{j} = importdata(file_vstm_alloc_5);
  %four_ibs_5{k} = importdata(file_four_ibs_5);
  %four_ibs_alloc_5{k} = importdata(file_four_alloc_5);
  %psy_ibs_5{k} = importdata(file_psy_ibs_5);
  %psy_ibs_alloc_5{k} = importdata(file_psy_alloc_5);
end
%%
vstm_nll_diff_arr = {};
vstm_nll_diff_mean_arr = {};
%four_nll_diff_arr = {};
%four_nll_diff_mean_arr = {};
psy_nll_diff_arr = {};
psy_nll_diff_mean_arr = {};

psy10_nll_diff_arr = {};
psy10_nll_diff_mean_arr = {};
for i = 1:120
    
    %four_nll_diff_cell = four_ibs_5{i} - four_ibs_5{i};
    %four_nll_diff_mean = mean(four_nll_diff_cell);
    %four_nll_diff_arr{end+1} = four_nll_diff_cell;
    %four_nll_diff_mean_arr{end+1} = four_nll_diff_mean

    psy_nll_diff_cell = psy_ibs_5{i} - psy_ibs_alloc_5{i}; %nll_exact difference 
    psy_nll_diff_mean = mean(psy_nll_diff_cell);
    psy_nll_diff_arr{end+1} = psy_nll_diff_cell;
    psy_nll_diff_mean_arr{end+1} = psy_nll_diff_mean
    
    psy10_nll_diff_cell = psy_ibs_10{i} - psy_ibs_alloc_10{i};
    psy10_nll_diff_mean = mean(psy10_nll_diff_cell);
    psy10_nll_diff_arr{end+1} = psy10_nll_diff_cell;
    psy10_nll_diff_mean_arr{end+1} = psy10_nll_diff_mean

end 
for i = 1:80
    vstm_nll_diff_cell = vstm_ibs_5{i}(1:35) - vstm_ibs_alloc_5{i}(1:35);
    vstm_nll_diff_mean = mean(vstm_nll_diff_cell);
    vstm_nll_diff_arr{end+1} = vstm_nll_diff_cell;
    vstm_nll_diff_mean_arr{end+1} = vstm_nll_diff_mean
end 

%%
vstm_result_mean_5 = mean(cell2mat(vstm_nll_diff_mean_arr))
psy_result_mean_5 = mean(cell2mat(psy_nll_diff_mean_arr))
psy_result_mean_10 = mean(cell2mat(psy10_nll_diff_mean_arr))

%%
output_psy_ibs_5 = {}
nll_psy_ibs_alloc5 = {}

output_vstm_ibs_5 = {}
nll_vstm_ibs_alloc5 = {}
for k = 1:numfiles
  %file_vstm_alloc_5 = sprintf('nll_exact_vstm_ibs_alloc_5_%d.txt', k);
  %file_four_alloc_5 = sprintf('nll_exact_fourinarow_ibs_alloc_5_%d.txt', k);
  output_psy_alloc_5 = sprintf('output_psycho_ibs_alloc_5_%d.txt', k);
  output_psy_ibs_alloc5{k} = importdata(output_psy_alloc_5);
  nll_psy_ibs_alloc5{k} =  output_psy_ibs_alloc5{k}(:,1);
  samples_psy_ibs_alloc5{k} =  output_psy_ibs_alloc5{k}(:,3);

  output_psy_5 = sprintf('output_psycho_ibs_5_%d.txt', k);
  output_psy_ibs_5{k} = importdata(output_psy_5);
  nll_psy_ibs_5{k} =  output_psy_ibs_5{k}(:,1);
  samples_psy_ibs_5{k} =  output_psy_ibs_5{k}(:,3);
end
for j = 1:80
  output_vstm_alloc_5 = sprintf('output_vstm_ibs_alloc_5_%d.txt', j);
  output_vstm_ibs_alloc5{j} = importdata(output_vstm_alloc_5);
  nll_vstm_ibs_alloc5{j} =  output_vstm_ibs_alloc5{j}(:,1);
  samples_vstm_ibs_alloc5{j} =  output_vstm_ibs_alloc5{j}(:,3);
  output_vstm_5 = sprintf('output_vstm_ibs_5_%d.txt', j);
  output_vstm_ibs_5{j} = importdata(output_vstm_5);
  nll_vstm_ibs_5{j} =  output_vstm_ibs_5{j}(:,1);
  samples_vstm_ibs_5{j} =  output_vstm_ibs_5{j}(:,3);


end
%% plot psycho nll results
%map = brewermap(3,'Set1'); 

mat_nll_psy_ibs_5 = cell2mat(nll_psy_ibs_5);
arr_nll_psy_ibs_5 = mat_nll_psy_ibs_5(:)';
mat_exact_psy_ibs_5 = cell2mat(psy_ibs_5);
arr_exact_psy_ibs_5 = mat_exact_psy_ibs_5(:)';

mat_nll_psy_ibs_alloc5 = cell2mat(nll_psy_ibs_alloc5);
arr_nll_psy_ibs_alloc5 = mat_nll_psy_ibs_alloc5(:)';
%mat_exact_psy_ibs_alloc5 = cell2mat(psy_ibs_alloc_5);
%arr_exact_psy_ibs_alloc5 = mat_exact_psy_ibs_alloc5(:)';



%xbar = [mat_nll_psy_ibs_5(1,:), mat_exact_psy_ibs_5(1,:)];
figure
histogram(arr_nll_psy_ibs_5,'edgecolor','none');
alpha(.3)
hold on
histogram(arr_exact_psy_ibs_5,'edgecolor','none');
alpha(.3)

hold on 
histogram(arr_nll_psy_ibs_alloc5,'edgecolor','none');

alpha(.3)
legend('nll ibs','nll exact','ibs alloc')
title('nll psycho ibs vs. exact vs alloc')
%iba alloc 5 vs. exact
%figure
%histogram(arr_nll_psy_ibs_alloc5,'edgecolor','none');
%alpha(.3)
%hold on
%histogram(arr_exact_psy_ibs_alloc5,'edgecolor','none');
%alpha(.3)
%legend('nll ibs alloc','nll exact')
%title('psycho ibs vs. exact')
%bar(mat_exact_psy_ibs_5(1,:));

%% vstm
rows = cellfun(@numel,nll_vstm_ibs_5);
cols = size(nll_vstm_ibs_5,2);
mat_nll_vstm_ibs_5 = nan(max(rows),cols);
for k = 1:cols
    mat_nll_vstm_ibs_5(1:rows(k),k) = nll_vstm_ibs_5{k};
end
%mat_nll_vstm_ibs_5 = cell2mat(nll_vstm_ibs_5);
arr_nll_vstm_ibs_5 = mat_nll_vstm_ibs_5(:)';

%exact
rows = cellfun(@numel,vstm_ibs_5);
cols = size(vstm_ibs_5,2);
mat_exact_vstm_ibs_5 = nan(max(rows),cols);

for k = 1:cols
    mat_exact_vstm_ibs_5(1:rows(k),k) = vstm_ibs_5{k};
end
%mat_exact_vstm_ibs_5 = cell2mat(vstm_ibs_5);
arr_exact_vstm_ibs_5 = mat_exact_vstm_ibs_5(:)';
%nll_alloc_5
rows = cellfun(@numel,nll_vstm_ibs_alloc5);
cols = size(nll_vstm_ibs_alloc5,2);
mat_nll_vstm_ibs_alloc5 = nan(max(rows),cols);
for k = 1:cols
    mat_nll_vstm_ibs_alloc5(1:rows(k),k) = nll_vstm_ibs_alloc5{k};
end
%mat_exact_vstm_ibs_5 = cell2mat(vstm_ibs_5);
arr_nll_vstm_ibs_alloc5 = mat_nll_vstm_ibs_alloc5(:)';
figure
histogram(arr_nll_vstm_ibs_5(1,:),'edgecolor','none');
alpha(.3)
hold on
histogram(arr_exact_vstm_ibs_5(1,:),'edgecolor','none');
alpha(.3)
hold on
histogram(arr_nll_vstm_ibs_alloc5,'edgecolor','none');
alpha(.3)
legend('nll ibs','nll exact','nll ibs alloc')
title('vstm ibs vs. exact vs. alloc')

%% samples used psycho
mat_samples_psy_ibs_5 = cell2mat(samples_psy_ibs_5);
arr_samples_psy_ibs_5 = mat_samples_psy_ibs_5(:)';
mat_samples_psy_ibs_alloc5 = cell2mat(samples_psy_ibs_alloc5);
arr_samples_psy_ibs_alloc5 = mat_samples_psy_ibs_alloc5(:)';

figure
histogram(arr_samples_psy_ibs_5,'edgecolor','none');
alpha(.3)
hold on 
histogram(arr_samples_psy_ibs_alloc5,'edgecolor','none');

alpha(.3)
legend('nll ibs','ibs alloc')
title('samples used psycho ibs vs. exact vs alloc')
%% samples used vstm
%mat_samples_vstm_ibs_5 = cell2mat(samples_vstm_ibs_5);
rows = cellfun(@numel,samples_vstm_ibs_5);
cols = size(samples_vstm_ibs_5,2);
mat_samples_vstm_ibs_5 = nan(max(rows),cols);

for k = 1:cols
    mat_samples_vstm_ibs_5(1:rows(k),k) = samples_vstm_ibs_5{k};
end
arr_samples_vstm_ibs_5 = mat_samples_vstm_ibs_5(:)';
%mat_samples_vstm_ibs_alloc5 = cell2mat(samples_vstm_ibs_alloc5);
rows = cellfun(@numel,samples_vstm_ibs_alloc5);
cols = size(samples_vstm_ibs_alloc5,2);
mat_samples_vstm_ibs_alloc5 = nan(max(rows),cols);

for k = 1:cols
    mat_samples_vstm_ibs_alloc5(1:rows(k),k) = samples_vstm_ibs_alloc5{k};
end
arr_samples_vstm_ibs_alloc5 = mat_samples_vstm_ibs_alloc5(:)';

figure
histogram(arr_samples_vstm_ibs_5,'edgecolor','none');
alpha(.3)
hold on 
histogram(arr_samples_vstm_ibs_alloc5,'edgecolor','none');

alpha(.3)
legend('nll ibs','ibs alloc')
title('samples used vstm ibs vs. exact vs alloc')