function all_in_one(c2)
%% encoding
m = 4;
q = 4;
h = log2(q);
cosNum = size(c2, 1);
w = log2(cosNum);
L = w + h*(m+1);
symNum = 100000;
data = randi([0, 1], symNum*L, 1);

c2 = gen_c2Matrix(c2);
[dataEncoded, G01, rep] = golay_encode(data, c2, q, w, m);
sk = exp(1j*2*pi/q * dataEncoded);
sk = reshape(sk, 2^m, []);
sn = ifft(sk);       

%% decoding
% golay_BER(data, sn, G01, rep, c2, m, h, L, {'mlReal', 'mlAbs', 'scatterReal', 'scatterAbs', 'bitwise'}, '2paths');
% golay_BER(data, sn, G01, rep, c2, m, h, L, {'mlReal'}, '2path');
Alamouti_BER(data, sk, G01, rep, c2, m, h, L,{'mlReal', 'mlAbs', 'scatterReal', 'scatterAbs', 'bitwise'});
if cosNum == 32
    Alamouti_mimo_BER(data, sk, G01, rep, c2, m, h, L, {'mlReal', 'mlAbs', 'scatterReal', 'scatterAbs', 'bitwise'}, 2);
    Alamouti_mimo_BER(data, sk, G01, rep, c2, m, h, L, {'mlReal', 'mlAbs', 'scatterReal', 'scatterAbs', 'bitwise'}, 4);
end
end