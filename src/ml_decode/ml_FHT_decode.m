function bits = ml_FHT_decode(SigRec, rep, q, type)
repNum = size(rep, 1);
repBitsLen = log2(repNum);
maxAbs = 0;
J = 1;
for j = 1:repNum
    SigRM1 = SigRec .* exp(-1j*2*pi/q*rep(j, :)).'; % additions equivalent for q = 2 or 4
    [M, bitsTemp] = RM1_qaryFHT_decode(SigRM1, q, type);
    if M > maxAbs
        maxAbs = M;
        RM1Bits = bitsTemp;
        J = j;
    end
end
repBits = de2bi(J-1, repBitsLen, 'left-msb').';
bits = [repBits; RM1Bits];
end
