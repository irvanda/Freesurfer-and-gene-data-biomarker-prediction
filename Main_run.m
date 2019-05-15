%%
clear; close all; clc; 
load data_using.mat;

%%
% Experiment 1: predicting freeSurfer biomarkers
t = cputime;
disp(['========================================']);
disp(['Predicting freeSurfer imaging markers starts']);
X_in = X(:, idx_freeSurf);
Y_in = Y_freeSurf(idx_freeSurf, :)';
group_idx = zeros(size(X, 1), 2);
group_idx(:, 1) = X_snp_gene_groups;
group_idx(:, 2) = group_info_r2{1}(:, 3);
l21Para = zeros(1, 2);
l21Para(1) = 1e-5; % para for group regularization
l21Para(2) = 1e-5; % para for L21-norm regularization
ridgePara = 1e-4;
[perf_out_freeSurf, B_our_freeSurf] = f_cal_lr_method_all_selected_top_k_multi_group(X_in, Y_in, group_idx, l21Para, ridgePara);
f_plot_performance(perf_out_freeSurf(1 : end - 1), 'FreeSurfer');
disp(['Predicting freeSurfer imaging markers is done by ', num2str(cputime - t, '%07.2f'), ' seconds.']);
set(gcf,'Position',[100 60 900 900],'PaperPositionMode','auto');
saveas(gcf,'Example01.tif','tif');

%%
% Experiment 2: predicting freeSurfer biomarkers
t = cputime;
disp(['========================================']);
disp(['Predicting VBM imaging markers starts']);
X_in = X(:, idx_vbm);
Y_in = Y_freeSurf(idx_vbm, :)';
group_idx = zeros(size(X, 1), 2);
group_idx(:, 1) = X_snp_gene_groups;
group_idx(:, 2) = group_info_r2{1}(:, 3);
l21Para = zeros(1, 2);
l21Para(1) = 1e-5; % para for group regularization
l21Para(2) = 1e-5; % para for L21-norm regularization
ridgePara = 1e-4;
[perf_out_vbm, B_our_vbm] = f_cal_lr_method_all_selected_top_k_multi_group(X_in, Y_in, group_idx, l21Para, ridgePara);
f_plot_performance(perf_out_vbm(1 : end - 1), 'VBM');
disp(['Predicting VBM imaging markers is done by ', num2str(cputime - t, '%07.2f'), ' seconds.']);
set(gcf,'Position',[100 60 900 900],'PaperPositionMode','auto');
saveas(gcf,'Example02.tif','tif');

%%
% Experiment 3: select SNP by using FreeSurfer imaging markers, group by gene
clear all; 
load data_using.mat;
X_in = X(:, idx_freeSurf);
Y_in = Y_freeSurf(idx_freeSurf, :)';
group_idx = X_snp_gene_groups;
l21Para = zeros(1, 2);
l21Para(1) = 1e-5; % para for group regularization
l21Para(2) = 2e-5; % para for L21-norm regularization
gene_names = cell(length(X_snp_gene_groups), 1);
for k = 1 : length(gene_names)
	gene_names{k} = gene_ids{X_snp_gene_groups(k)};
end
[snp_info, B_raw] = f_cal_snp_selection(X_in, Y_in, group_idx, snp_id, gene_ids, l21Para, gene_names, 1);
clearvars -except snp_info;

[x,y]=size(snp_info);
fout=fopen('Example03.csv','w');
for i=1:x
    for j=1:y
        if strcmpi(class(snp_info{i,j}),'char')
            fprintf(fout,'%s,',snp_info{i,j});
        elseif strcmpi(class(snp_info{i,j}),'double')
            fprintf(fout,'%.3e,',snp_info{i,j});
        end
    end
    fprintf(fout,'\n');
end
fclose(fout);