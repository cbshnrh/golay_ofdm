function bits = bitwise_decode(r, m, h, rep, G01)
repNum = size(rep, 1);
idxEncoding = 1:repNum;
k = 0;
W = zeros(1, m+1);
f = zeros(2^m, h);
while k <= h-1
    dist = distinct_rep(rep, 2^(k+1));
    z = dist.rep;
    s = size(z, 1);
    p = 1;
    Y = 0;
    while p <= s
        yn = 2^(k-1) - weight_specific(r - z(p, :).', 2^(k+1));
        yh = fwht(yn, [], 'hadamard');
        [v, j] = max(abs(yh));
        if v > abs(Y)
            Y = yh(j);
            J = j-1;
            P = p;
        end
        p = p+1;
    end
    if Y > 0
        w0 = 0;
    else
        w0 = 1;
    end
    w = de2bi(J, m, 'left-msb');
    W = W + 2^k*[w0, w];
    f(:, k+1) = mod([w0, w]*G01, 2^(h-k)).';
    idx = dist.idx;
    rep = rep(idx == P, :);
    idxEncoding = idxEncoding(idx == P);
    r = mod(r-2^k*f(:, k+1), 2^h);
    k = k+1;
end
bitsRep = de2bi(idxEncoding-1, log2(repNum), 'left-msb');
bitsRM1 = de2bi(W, h, 'left-msb').';
bits = [bitsRep.'; bitsRM1(:)];
end

function dist = distinct_rep(rep, m)
repMod = mod(rep, m);
[c, ~, ic] = unique(repMod, 'rows', 'stable');
dist.rep = c;
dist.idx = ic;
end

function w = weight_specific(r, m)
w = min(mod(r, m), m-mod(r, m));
end