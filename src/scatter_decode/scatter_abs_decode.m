function bits = scatter_abs_decode(r, c2, rep, q)
m = log2(length(r));
rk = gen_rk(r);
L = size(c2, 3);
index = gen_index(c2, q);

    
Y = zeros((q/2)^(m-1), m);
for k = 1:m
    Y(:, k) = qary_fht(rk(:, k), q/2);
end

d = zeros(L, 1);
for i = 1:L
    for k = 1:m
        d(i) = d(i) + abs(Y(index(i, k), k))^2;
    end
end
[~, I] = max(d);

c01 = zeros(m+1, 1);
for k = 1:m
    RP = 0;
    for j = 0:q-1
        realPart = real(exp(-1i*j*2*pi/q) * Y(index(I, k), k));
        if RP < realPart
            RP = realPart;
            ak = j;
        end
    end
    c01(k+1) = ak;
end

z = RM2_encode(rep(I, :), c01, q);
R = exp(-1j*2*pi*z/q) * r;

RP = 0;
for i = 0:q-1
    realPart = real(exp(-1j*i*2*pi/q) * R);
    if RP < realPart
        RP = realPart;
        c01(1) = i;
    end
end

repBits = de2bi(I-1, log2(L), 'left-msb').';
RM1Bits = de2bi(c01, log2(q), 'left-msb');
RM1Bits = reshape(RM1Bits.', [], 1);
bits = [repBits; RM1Bits];
end


            
            
        
