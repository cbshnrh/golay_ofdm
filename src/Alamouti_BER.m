function Alamouti_BER(data, sk, G01, rep, c2, m, h, L, types)
close all;
figure; 
flag = {'r-p', 'g-o', 'b-*', 'c-x', 'm-s', 'y-d'};
addpath(genpath(pwd));

PG = db2pow([-2.3 -3.8]);
symNum = size(sk, 2);

for i = 1:length(types)
    type = types{i};
    init = 8;
    EbN0 = init;
    ber = 1;
    E = [];
    step = 1;
    while ber >= 1e-4
        hCoff1 = sqrt(1/2) * repmat(sqrt(PG), symNum/2, 1)...
            .* (randn(symNum/2, 2) + 1j*randn(symNum/2, 2));
        hCoff2 = sqrt(1/2) * repmat(sqrt(PG), symNum/2, 1)...
            .* (randn(symNum/2, 2) + 1j*randn(symNum/2, 2));
        EbN0
        SNR = db2pow(EbN0-3-10*log10((2^m)/L));
        [sk1, sk2] = gen_diversity(sk);
        sn1 = ifft(sk1);
        sn2 = ifft(sk2);
%         PAPR(sn1(:, 1))
%         PAPR(sn2(:, 1))
        [sn1Fading, lambda1] = ofdm_multipath(sn1, hCoff1);
        [sn2Fading, lambda2] = ofdm_multipath(sn2, hCoff2);
        snFading = sn1Fading + sn2Fading;
        snPow = mean(abs(snFading).^2, 'all');
        noise = sqrt(snPow/SNR) * 1/sqrt(2) * (randn(size(snFading)) + 1j*randn(size(snFading)));
        snRec = snFading + noise;
        skRec = fft(snRec);
        skRec(:, 2:2:end) = conj(skRec(:, 2:2:end));
        % solve linear functions
        for j = 1:size(skRec, 1)
            for k = 1:2:size(skRec, 2)
                h1 = lambda1(j, (k-1)/2+1);
                h2 = lambda2(j, (k-1)/2+1);
                H = [h1, h2; conj(h2), -conj(h1)];
                skRec(j, k:k+1) = skRec(j, k:k+1) * conj(H);
            end
        end
        skRec = skRec(:);

        if isequal('bitwise', type)
            skRec = angle(skRec)*2^h/(2*pi);
            skRec = mod(round(skRec), 2^h);
        end
        bits = golay_decode(skRec, G01, rep, c2, m, h, L, type, EbN0);
        ber = length(find(data - bits)) / length(data)
        E = [E, ber];
        EbN0 = EbN0 + step;
    end
    semilogy(init:step:EbN0-step, E, flag{i});
    hold on;   
    grid on;
end
legend(types);
ylabel('BER');
xlabel('E_b/N_0 (dB)');
title(sprintf('m = 4, %d cosets, PG = [-2.3 -3.8], Alamouti, QPSK', size(rep, 1)));
savefig(sprintf('./figure/2paths_%dcosets_alamouti', size(rep, 1)));
end

function [sk1, sk2] = gen_diversity(sk)
sk1 = sk;
sk2 = zeros(size(sk));
for i = 1:2:size(sk1, 2)
    sk1(:, i+1) = -conj(sk1(:, i+1));
    sk2(:, i) = sk(:, i+1);
    sk2(:, i+1) = conj(sk(:, i));
end
end

function [snFading, lambda] = ofdm_multipath(sn, hCoff)
[N, M] = size(sn);
snFading = zeros(size(sn));
for i = 1:M/2
    hn = [hCoff(i, 1), hCoff(i, 2), zeros(1,N-2)];
    snFading(:, 2*i-1) = cconv(hn.', sn(:, 2*i-1), N);
    snFading(:, 2*i) = cconv(hn.', sn(:, 2*i), N);
end
lambda = fft(cat(1, hCoff.', zeros(N-2, M/2)));
end

function R = PAPR(sn)
envPower = sn .* conj(sn);
peak = max(envPower);
avg = mean(envPower);
R = peak/avg;
end