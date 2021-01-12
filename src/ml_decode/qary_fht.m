function Y = qary_fht(y, q)
m = log2(length(y));
xi = exp(1j*2*pi/q);
H1 = ones(2, q);
H1(2, :) = conj(xi).^(0:q-1);
Y = y;
for i = 1:m
    HSparse = kron(eye(2^(m-i)), H1);
    HSparse = kron(HSparse, eye(q^(i-1)));
    Y = HSparse.' * Y;
end
end