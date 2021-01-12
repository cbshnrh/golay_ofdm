function bits= golay_decode(snRec, G01, rep, c2, m, h, L, type, EbN0)
N = length(snRec) / 2^m;
bits = zeros(L*N, 1);

wb = waitbar(0, sprintf('Eb/N0 = %0.1f dB', EbN0));
for i = 1:N
    waitbar(i/N, wb)
    r = snRec(1+(i-1)*2^m : i*2^m);
    if isequal(type, 'bitwise')
        bits(1+(i-1)*L : i*L) = bitwise_decode(r, m, h, rep, G01); 
    elseif isequal(type, 'mlAbs')
        bits(1+(i-1)*L : i*L) = ml_FHT_decode(r, rep, 2^h, 'abs');
    elseif isequal(type, 'mlReal')
        bits(1+(i-1)*L : i*L) = ml_FHT_decode(r, rep, 2^h, 'real');
    elseif isequal(type, 'scatterAbs')
        bits(1+(i-1)*L : i*L) = scatter_abs_decode(r, c2, rep, 2^h);
    elseif isequal(type, 'scatterReal')
        bits(1+(i-1)*L : i*L) = scatter_real_decode(r, c2, rep, 2^h);
    end
end
delete(wb);
end
