function golay_BER(data, sn, G01, rep, c2, m, h, L, types, channel)
close all;  
figure;
flag = {'r-p', 'g-o', 'b-*', 'c-x', 'm-s', 'y-d'};
addpath(genpath(pwd));
PG = db2pow([-2.3 -3.8]);
symNum = size(sn, 2);


for i = 1:length(types)
    type = types{i};
    init = 8;
    EbN0 = init;
    ber = 1;
    E = [];
    step = 1;
    while ber >= 1e-4
        EbN0
        if isequal(channel, '2paths')
            hCoff = sqrt(1/2) * repmat(sqrt(PG), symNum, 1)...
                 .* (randn(symNum, 2) + 1j*randn(symNum, 2));
            snFading = zeros(2^m, symNum);
            lambda = zeros(2^m, symNum);
            for j = 1:symNum
                hn = [hCoff(j, 1), hCoff(j, 2), zeros(1, 2^m-2)];
                snFading(:, j) = cconv(hn.', sn(:, j), 2^m);
                lambda(:, j) = fft([hCoff(j, 1), hCoff(j, 2), zeros(1, 2^m-2)]).';
            end
            snPow = mean(abs(snFading).^2, 'all');
            SNR = db2pow(EbN0-10*log10(2^m/L));
            noise = sqrt(snPow/SNR) * 1/sqrt(2) * (randn(size(snFading)) + 1j*randn(size(snFading)));
            snRec = snFading + noise;
            skRec = fft(snRec);
            skRec = skRec .* conj(lambda);
        else
            snRec = awgn(sn, EbN0-10*log10(2^m/L), 'measured');
            skRec= fft(snRec);
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
if isequal(channel, 'awgn')
    tit = sprintf('m = 4, %d cosets, awgn, QPSK', size(rep, 1));
    filename = sprintf('./figure/awgn_%dcosets', size(rep, 1));
else
    tit = sprintf('m = 4, %d cosets, PG = [-2.3 -3.8], QPSK', size(rep, 1));
    filename = sprintf('./figure/2paths_%dcosets', size(rep, 1));
end
title(tit);
savefig(filename);