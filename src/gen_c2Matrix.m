function c2Matrix = gen_c2Matrix(c2)
[L, M] = size(c2);
syms x
solution = solve(x*(x-1)/2 == M, x);
M = eval(solution(2));
c2Matrix = zeros(M, M, L);
for i = 1:L
    index = 0;
    for j = 1:M-1
        for k = j+1:M
            index = index + 1;
            c2Matrix(j, k, i) = c2(i, index);
        end
    end
end
end