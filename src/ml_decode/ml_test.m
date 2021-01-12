%% qary_fht
q = 4;
m = 4;
G = de2bi(0:2^m-1, m, 'left-msb');
G = G.';
G = [ones(1, 2^m); G];
L = log2(q) * (m+1);
EbN0 = 5.6;
cycle = 10000;
e = 0;
for i = 1:cycle
    info = randi([0 1], m+1, log2(q));
    phase = bi2de(info, 'left-msb');
    y = mod(phase.' * G, q).';
    y = exp(1j*2*pi/q*y);
    y = awgn(y, EbN0-10*log10(2^m/L), 'measured');
    Y = qary_fht(y, q);
%     Y1 = qary_ht(y, q);
%     Y - Y1    

% RM1_qaryFHT_decode
    [M, bits] = RM1_qaryFHT_decode(y, q);
    info = reshape(info.', [], 1);
    errorRate = sum(bits ~= info) / length(bits);
    e = e + errorRate;
end
e = e/cycle