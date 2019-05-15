function [snp_info, B_raw] = f_cal_snp_selection(X, Y, group_idx, snp_name, group_name, l21Para, gene_name, gene_grouping)

% Input:
% 	X: p x n
% 	Y: c x n
% 	l21para: 2 x 1, the first is for group l21 norm, the second is for l21 norm
% Output:
% 	B: p x c

if size(X, 2) ~= size(Y, 2)
	error(dbstack);
end
if length(group_name) ~= length(unique(group_idx))
	error(dbstack);
end

if nargin < 5
	error(dbstack);
end
if nargin < 6
	l21Para = [10, 10];
end

x_normalized = f_cal_normalized_feature(X);
y_normalized = f_cal_normalized_feature(Y);

t = cputime;
%disp(['Our method computation starts ...']);
B_raw = L2R21R21group(x_normalized', y_normalized', group_idx, l21Para(1), l21Para(2));
%disp(['Our method computation is done by ', num2str(cputime - t, '%06.2f'), ' seconds']);

%[snp_info, group_info] = f_process_B(B_raw, group_idx, snp_name, group_name);
[snp_info] = f_process_B_new1(B_raw, group_idx, snp_name, gene_name, gene_grouping);

