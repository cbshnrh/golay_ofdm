function rep = gen_rep(c2, G01, q)
[M, M, L] = size(c2);
rep = zeros(L, 2^M);
for i = 1:L
    [I, J] = ind2sub([M, M], find(c2(:, :, i)));
    for k= 1:length(I)
        rep(i, :) = rep(i, :) + 2 * G01(I(k)+1, :) .* G01(J(k)+1, :);
    end
    rep(i, :) = mod(rep(i, :), q);
end
end
