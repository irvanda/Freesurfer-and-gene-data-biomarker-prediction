function [group_idx, num_groups, group_stat] = f_cal_construct_group_frm_w(W, thres)

if ~isequal(W, W')
	error('Input w is asymmetric');
end

if nargin < 2
	thres = 0.4;
end

num_thres = length(thres);

group_idx = zeros(size(W, 1), num_thres);
num_groups = zeros(num_thres, 1);
group_stat = cell(num_thres, 1);
for k_thres = 1 : num_thres
	W_using = W;
	W_using(W_using < thres(k_thres)) = 0;
	[num_groups(k_thres), group_idx(:, k_thres)] = graphconncomp(W_using);

	group_stat{k_thres} = zeros(num_groups(k_thres), 2);
	for k = 1 : num_groups(k_thres)
		group_stat{k_thres}(k, 1) = k;
		group_stat{k_thres}(k, 2) = nnz(group_idx(:, k_thres) == k);
	end
end

group_stat{num_thres + 1} = zeros(num_thres * 2);
