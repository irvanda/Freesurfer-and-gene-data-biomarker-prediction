function [X_normalized] = f_cal_normalized_feature(X, isRow)

% We do L1 normalization and centralization.
% X: p x n

if nargin < 2
	isRow = 1;
end

if isRow
	X_normalized = diag(1 ./ sum(abs(X), 2)) * X;
	X_normalized = X_normalized - repmat(mean(X_normalized, 2), 1, size(X_normalized, 2));
else
	X_normalized = X * diag(1 ./ sum(abs(X), 1));
	X_normalized = X_normalized - repmat(mean(X_normalized, 1), size(X_normalized, 1), 1);
end


