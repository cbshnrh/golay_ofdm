%% encoding
m = 4;
q = 4;
h = log2(q);
cosNum = 32;
w = log2(cosNum);
L = w + h*(m+1);
symNum = 10000;
data = randi([0, 1], symNum*L, 1);
c2 = [
    0 0 1 1 0 1; 0 1 0 0 1 1; 0 1 0 1 1 0; 0 1 1 1 0 0;
    1 0 0 0 1 1; 1 0 0 1 0 1; 1 0 1 0 0 1; 1 0 1 1 0 0;
    1 1 0 0 0 1; 1 1 0 0 1 0; 0 0 1 1 1 0; 0 1 1 0 1 0;
    0 0 0 1 0 1; 0 0 0 1 1 0; 0 0 1 0 1 1; 0 0 1 1 0 0;
    0 0 1 1 1 1; 0 1 0 0 0 1; 0 1 0 0 1 0; 0 1 0 1 0 1;
    0 1 0 1 1 1; 0 1 1 0 0 0; 0 1 1 0 1 1; 0 1 1 1 0 1;
    0 1 1 1 1 1; 1 0 0 0 0 1; 1 0 0 0 1 0; 1 0 0 1 0 0;
    1 0 0 1 1 0; 1 0 0 1 1 1; 1 0 1 0 0 0; 1 0 1 0 1 1];
c2 = gen_c2Matrix(c2);
[dataEncoded, G01, rep] = golay_encode(data, c2, q, w, m);
sk = exp(1j*2*pi/q * dataEncoded);
sk = reshape(sk, 2^m, []);
sn = ifft(sk);
sn = cat(1, sn, zeros(1, size(sn, 2)));
sn = sn(:);


%% decoding
errorCorrectSim(data, sn, G01, rep, c2, m, h, L, {'bitwise', 'scatterAbs', 'scatterReal', 'mlAbs', 'mlReal'}, '2path');
% errorCorrectSim(data, sn, G01, rep, c2, m, h, L, {'scatterAbs', 'scatterReal'});
% errorCorrectSim(data, sn, G01, rep, c2, m, h, L, {'bitwise'});
