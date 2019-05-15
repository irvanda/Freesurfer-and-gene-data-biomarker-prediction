function [snp_info] = f_process_B(B, group_idx, snp_name, gene_name, grouping_by_gene)

if nargin < 4
	error('Insufficient para.');
end
if nargin < 5
	grouping_by_gene = 1;
end

B_weight = sum(abs(B), 2);
[B_weight_sorted, B_weight_sorted_idx] = sort(B_weight, 'descend');

group_ids = unique(group_idx);
for k = 1 : length(group_ids)
	group_weight(k) = sum(B_weight(group_idx == group_ids(k)));
end

snp_info = cell(size(B, 2) + 5, size(B, 1));
for k = 1 : size(B, 1)
	snp_info{1, k} = snp_name{B_weight_sorted_idx(k)};
	snp_info{2, k} = B_weight(B_weight_sorted_idx(k));
	if grouping_by_gene
		snp_info{3, k} = gene_name{B_weight_sorted_idx(k)};
	else
		snp_info{3, k} = group_idx(B_weight_sorted_idx(k));
	end
	snp_info{4, k} = group_weight(group_idx(B_weight_sorted_idx(k)));
	snp_info{5, k} = [snp_info{1, k}, '-', gene_name{B_weight_sorted_idx(k)}];
	for kk = 1 : size(B, 2)
		snp_info{5 + kk, k} = B(B_weight_sorted_idx(k), kk);
	end
end


