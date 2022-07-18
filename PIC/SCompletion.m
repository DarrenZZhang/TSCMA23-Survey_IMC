function S = SCompletion(W, index, issymmetric)
if nargin < 3
    issymmetric = 1;
end;

n_view = length(W);
S = cell(1, n_view);
[N, ~] = size(W{1});
for v = 1:n_view
    S{v} = zeros(N, N);
    for id = 1:N
        v_count = false;
        if (find(index{v} == id))
            S{v}(id, :) = W{v}(id, :);
        else          
            for vv = 1:n_view
                if (find(index{vv} == id))
                    S{v}(id, :) = S{v}(id, :) + W{vv}(id, :);
                    v_count = true;
                end
            end
        end
        if v_count
            S{v}(id, :) = S{v}(id, :)/sum(S{v}(id, :));
        end
    end
end
if 1 == issymmetric
    for v = 1:n_view
        S{v} = (S{v}+S{v}')/2;
    end
end
end