function bits = scatter_real_decode(r, c2, rep, q)
m = log2(length(r));
rk = gen_rk(r);
L = size(c2, 3);
index = gen_index(c2, q);

    
Y = zeros((q/2)^(m-1), m);
for k = 1:m
    Y(:, k) = qary_fht(rk(:, k), q/2);
end
Y3 = cat(3, Y, zeros([size(Y), q-1]));
for i = 1:q-1
    Y3(:, :, i+1) = Y * exp(-1j*2*pi/q * i);
end
[YRealMax, Ind] = max(real(Y3), [], 3);

d = zeros(L, 1);
for i = 1:L
    for k = 1:m
        d(i) = d(i) + YRealMax(index(i, k), k);
    end
end

[~, I] = max(d);
c01 = zeros(m+1, 1);
for k = 1:m
    c01(k+1) = Ind(index(I, k), k) - 1;
end

z = RM2_encode(rep(I, :), c01, q);
R = exp(-1j*2*pi*z/q) * r;
R = R * exp(-1j*2*pi/q * (0:q-1));
[~, c0] = max(real(R));
c01(1) = c0 - 1;

repBits = de2bi(I-1, log2(L), 'left-msb').';
RM1Bits = de2bi(c01, log2(q), 'left-msb');
RM1Bits = reshape(RM1Bits.', [], 1);
bits = [repBits; RM1Bits];
end


            
            
        
