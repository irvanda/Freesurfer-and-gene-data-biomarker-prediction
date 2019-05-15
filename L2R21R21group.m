function [X, obj]=L2R21R21group(A, Y, class_idx, r1, r2)
%% 2-norm loss with general 21-norm regularization

%% Problem
%
%  min_X  || A X - Y||^2 + r1 * ||X||_21group + r2*||X||_21



[dim, n] = size(A);
class_set = unique(class_idx);
class_num = length(class_set);

AA = A'*A;
Ay = A'*Y;
d1 = ones(n,1);
d2 = ones(n,1);
Xi = zeros(n, 1);

for iter = 1:100
    D1 = diag(d1);
    D2 = diag(d2);
    X = (AA+r1*D1+r2*D2)\Ay;
    
    for c = 1:class_num
        idx = find(class_idx==class_set(c));
        Xc = X(idx,:);
        di = sqrt(sum(sum(Xc.*Xc))+eps);
        Xi(idx) = di;
        ob(c) = di;
    end;
    d1 = 0.5./Xi;
    
    Xi = sqrt(sum(X.*X,2)+eps);
    d2 = 0.5./(Xi);
    ob2 = sum(Xi);
    
    obj(iter) = trace((A*X-Y)'*(A*X-Y)) + r1*sum(ob) + r2*ob2;
end;

