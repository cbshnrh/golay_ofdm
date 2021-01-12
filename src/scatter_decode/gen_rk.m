function rk = gen_rk(r)
m = log2(length(r));
rk = zeros(2^(m-1), m);
IK = gen_IK(m);
for k = 1:m
    rk(:, k) = r(IK(k, :)+1 + 2^(m-k)) .* conj(r(IK(k, :)+1));
end
end