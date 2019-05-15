function f_plot_performance(perf_out, title_txt, top_k_selected, legend_txt)

if nargin < 2
	title_txt = 'FreeSurfer';
end
if nargin < 3
	top_k_selected = [10 : 10 : 100];
end
if nargin < 4
	legend_txt = {
		'Multivariate regression';
		'Ridge regression';
		'MTFL';
		'Our method (gene)';
		'Our method (r^2>0.2)';
	};
end

if length(perf_out) ~= length(top_k_selected)
	error(dbstack);
end

line_width = 2;
marker_size = 18;
font_size = 30;
line_style = {'-.bo', '-.rs', '-kx', '-bv', '--rp', '-b^', '--r*'};


num_points = length(perf_out);

rmse = zeros(num_points, size(perf_out{1}, 2) / 3);
for k = 1 : num_points
	for kk = 1 : size(perf_out{k}, 2) / 3
		rmse(k, kk) = perf_out{k}(end, 1 + (kk - 1) * 3);
	end
end

figure; hold on; box on; grid on;
for k = 1 : size(rmse, 2)
	plot(top_k_selected, rmse(:, k), line_style{k}, 'LineWidth', line_width, 'MarkerSize', marker_size);
end
legend(legend_txt, 'FontSize', font_size - 6, 'Location', 'NorthWest');
xlabel('Number of selected SNPs', 'FontSize', font_size);
ylabel('Rooted Mean Square Error', 'FontSize', font_size);
title([title_txt, ' imaging phenotypes'], 'FontSize', font_size + 4);
%xlim([10, 100]);
set(gca, 'FontSize', font_size - 4, 'LineWidth', line_width);
