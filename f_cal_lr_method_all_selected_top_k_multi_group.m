function [perf_out, B_our, data_stat] = f_cal_lr_method_all_selected_top_k_multi_group(X, Y, group_idx, l21Para, ridgePara, top_k_selected, idx_test)

% Input:
% 	X: p x n
% 	Y: c x n
% 	group_idx: p x k, where k the number of grouping measurements
% 	l21para: 2 x 1, the first is for group l21 norm, the second is for l21 norm
% Output:
% 	B: c x p

if size(X, 2) ~= size(Y, 2)
	error(dbstack);
end

if nargin < 4
	l21Para = [1e-5, 1e-5];
end
if nargin < 5
	ridgePara = 1e-4;
end
if nargin < 6
	top_k_selected = [10 : 10 : 100];
end
if nargin < 7
	% by default, we use 80% data for test
	indices = crossvalind('Kfold', size(X, 2), 5);
	idx_test = indices == 1;
end

x_normalized = f_cal_normalized_feature(X);
y_normalized = f_cal_normalized_feature(Y);

data_stat{1} = idx_test;
data_stat{2} = mean(x_normalized, 2);
data_stat{3} = std(x_normalized, 0, 2);

x_train = x_normalized(:, ~idx_test);
x_test = x_normalized(:, idx_test);
y_train = y_normalized(:, ~idx_test);
y_test = y_normalized(:, idx_test);

t = cputime;


% our method
disp(['Our method computation starts ...']);
for kk = 1 : size(group_idx, 2)
	B_our{kk} = L2R21R21group(x_train', y_train', group_idx(:, kk), l21Para(1), l21Para(2));
	idx_selected = f_select_feauture(B_our{kk}, top_k_selected);
	for k = 1 : length(top_k_selected)
		y_test_predicted_our{kk, k} = B_our{kk}(idx_selected(1 : top_k_selected(k)), :)' * x_test(idx_selected(1 : top_k_selected(k)), :);
		perf_our{kk, k} = f_cal_lr_perf((y_test_predicted_our{kk, k})', y_test');
	end
end
disp(['Our method computation (', num2str(length(idx_selected)), ' selected features) is done by ', num2str(cputime - t, '%06.2f'), ' seconds']);

% our method no group
disp(['Our method no group computation starts ...']);
B_our_nogroup = L2R21R21group(x_train', y_train', group_idx(:, 1), 0, 5e-7);
idx_selected = f_select_feauture(B_our_nogroup, top_k_selected);
for k = 1 : length(top_k_selected)
	y_test_predicted_our_nogroup{k} = B_our_nogroup(idx_selected(1 : top_k_selected(k)), :)' * x_test(idx_selected(1 : top_k_selected(k)), :);
	perf_our_nogroup{k} = f_cal_lr_perf((y_test_predicted_our_nogroup{k})', y_test');
end
disp(['Our method no group computation (', num2str(length(idx_selected)), ' selected features) is done by ', num2str(cputime - t, '%06.2f'), ' seconds']);

% linear regression
disp(['Linear regression computation starts ...']);
%B_lr = y_train * x_train' / (x_train * x_train');
B_lr =  y_train * x_train' * pinv(x_train * x_train');
idx_selected = f_select_feauture(B_lr', top_k_selected);
for k = 1 : length(top_k_selected)
	y_test_predicted_lr{k} = B_lr(:, idx_selected(1 : top_k_selected(k))) * x_test(idx_selected(1 : top_k_selected(k)), :);
	perf_lr{k} = f_cal_lr_perf((y_test_predicted_lr{k})', y_test');
end
disp(['Linear regression computation (', num2str(length(idx_selected)), ' selected features) is done by ', num2str(cputime - t, '%06.2f'), ' seconds']);

% ridge regression
disp(['Rigde regression computation starts ...']);
B_ridge = y_train * x_train' / (x_train * x_train' + ridgePara * eye(size(x_train, 1)));
idx_selected = f_select_feauture(B_ridge', top_k_selected);
for k = 1 : length(top_k_selected)
	y_test_predicted_ridge{k} = B_ridge(:, idx_selected(1 : top_k_selected(k))) * x_test(idx_selected(1 : top_k_selected(k)), :);
	perf_ridge{k} = f_cal_lr_perf((y_test_predicted_ridge{k})', y_test');
end

perf_out = cell(length(top_k_selected) + 1, 1);
perf_out{length(top_k_selected) + 1} = []; 
for k = 1 : length(top_k_selected)
	perf_out{k} = horzcat(perf_lr{k}, perf_ridge{k});
	perf_out{k} = horzcat(perf_out{k}, perf_our_nogroup{k});
	for kk = 1 : size(group_idx, 2)
		perf_out{k} = horzcat(perf_out{k}, perf_our{kk, k});
	end
	perf_out{length(top_k_selected) + 1} = vertcat(perf_out{length(top_k_selected) + 1}, perf_out{k});
end

%%%% end of f_cal_lr_method_all_selected %%%%



function idx_selected = f_select_feauture(B, top_k_selected)

% B: p x c

if nargin < 2 
	top_k_selected = 10;
end

b_weight = sum(abs(B), 2);
[b_weigth_sorted, b_weigth_sorted_idx] = sort(b_weight, 'descend');

idx_selected = b_weigth_sorted_idx(1 : max(top_k_selected));

%%%% end of f_select_feauture %%%%
