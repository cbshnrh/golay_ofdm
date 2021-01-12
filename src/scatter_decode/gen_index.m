function index = gen_index(c2, q)
[M, M, L] = size(c2);
index = zeros(L, M);
for i = 1:L
    for j = 1:M
        vec = [c2(1:j-1, j, i); c2(j, j+1:M, i).'];
        index(i, j) = (q/2).^(M-2:-1:0) * vec+ 1;
    end
end
end