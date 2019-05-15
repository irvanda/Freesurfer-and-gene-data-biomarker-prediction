function [perf_out] = f_cal_lr_perf(y_predicted, y_ground_truth)

% input
% 	y_predicted: n x c
% 	y_ground_truth: n x c
% output
% 	perf_out:
% 		1: rmse
% 		2~3: to be used

if ~isequal(size(y_predicted), size(y_ground_truth))
	error(dbstack);
end

[num_data, num_class] = size(y_predicted);

perf_out = zeros(num_class, 3);

for k = 1 : num_class
	perf_out(k, 1) = sqrt(norm(y_predicted(:, k) - y_ground_truth(:, k)) .^ 2 / num_data);
end

perf_out(num_class + 1, :) = mean(perf_out, 1);
