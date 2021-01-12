function z = RM2_encode(rep, c01, q)
m = length(c01) - 1;
G1 = de2bi(0:2^m-1, m, 'left-msb');
G1 = G1.';
z = mod(c01(2:end).' * G1 + rep, q);
end