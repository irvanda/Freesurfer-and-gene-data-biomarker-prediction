%clc
%clear
% Experiment 1: predicting freeSurfer biomarkers
t = cputime;
disp(['========================================']);
disp(['Predicting freeSurfer thickness imaging markers starts']);
X_in = randi([0 2],1224, 50);
Y_in = load('fs_thickness.txt');
ridgePara = 1e-4;
[perf_out_freeSurf] = f_cal_lr_method_all_selected_top_k_multi_group_try(X_in, Y_in, ridgePara);
f_plot_performance(perf_out_freeSurf(1 : end - 1), 'FreeSurfer_thickness');
disp(['Predicting freeSurfer imaging markers is done by ', num2str(cputime - t, '%07.2f'), ' seconds.']);
set(gcf,'Position',[100 60 900 900],'PaperPositionMode','auto');
saveas(gcf,'Example01.tif','tif');