function I = gen_IK(m)
I = zeros(m, 2^(m-1));
num = 0:2^m-1;
bits = de2bi(num, m, 'left-msb');
for i = 1:m
    I(i, :) = bi2de(bits(~bits(:, i), :), 'left-msb'); %% logic index
end
end
