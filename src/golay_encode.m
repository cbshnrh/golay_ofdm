function [dataEncoded, G01, reps] = golay_encode(data, c2, q, w, m)
h = log2(q);
L = w + h*(m+1);
N = ceil(length(data)/L);
num0 = N*L - length(data);
data = [data; zeros(num0, 1)];
G0 = de2bi(0:2^m-1, m, 'left-msb');
G0 = G0.';
G01 = [ones(1, 2^m); G0];
reps = gen_rep(c2, G01, q);
dataEncoded = zeros(N*2^m, 1);
for i = 1:N
    RM1Bits = data((w+1:L)+(i-1)*L);
    RM1Bits = reshape(RM1Bits, h, m+1);
    c01 = bi2de(RM1Bits.', 'left-msb');
    repBits = data((1:w)+(i-1)*L);
    repIndex = bi2de(repBits.', 'left-msb') + 1;
    rep = reps(repIndex, :);
    dataEncoded((1:2^m)+(i-1)*2^m) = mod(c01.'*G01 + rep, q);
end
end

    
