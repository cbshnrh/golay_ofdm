function [M, bits] = RM1_qaryFHT_decode(y, q,type)
h = log2(q);
m = log2(length(y));
assert(mod(h, 1) == 0, "q must be a power of 2");
Y = qary_fht(y, q);
if isequal(type, 'real')
    Y = Y * exp(-1j*2*pi/q * (0:q-1));
    [M2, I2] = max(real(Y), [], 2);
    [M, I] = max(M2);
    J = I2(I) - 1;
else
    [M, I] = max(abs(Y));
    phi = mod(angle(Y(I)), 2*pi);
    if pi*(2*q-1)/q < phi && phi <= 2*pi
        phi = phi - 2*pi;
    end

    for j = 0:q-1
        if pi*(2*j-1)/q < phi && phi <= pi*(2*j+1)/q
            J = j;
            break
        end
    end
end

c1 = de2bi(I-1, h*m, 'left-msb');
c0 = de2bi(J, h, 'left-msb');
bits = [c0.'; c1.'];
end