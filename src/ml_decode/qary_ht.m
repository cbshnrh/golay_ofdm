function Y = qary_ht(y, q)
m = log2(length(y));
Y = zeros(q^m, 1);
G = de2bi(0:2^m-1, m, 'left-msb');
G = G.';
for i = 0:q^m-1
    u = de2bi(i, log2(q^m), 'left-msb');
    u = reshape(u, log2(q), []).';
    u = bi2de(u, 'left-msb');
    ux = mod(u.' * G, q);
    Y(i+1) = exp(-1j*2*pi/q*ux) * y;
end
end