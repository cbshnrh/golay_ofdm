function ht = gen_channel(pathDelays,averagePathGains,sampleRate,winLen)
% pathDelays:多径的时延
% averagePathGains:多径的平均增益
% sampleRate:采样率  此处sampleRate=1(不用管)
% winLen:sinc函数的窗口大小
% ht:生成的时域多径信道(整数倍采样)
% winLen = 40;
averagePathPow = db2pow(averagePathGains);
taud = pathDelays*sampleRate;
rollOffFactor = 0.8;
upSampRate = 1;
h = zeros(1,winLen*2*upSampRate+1);
% genieChan = zeros(1,length(pathDelays));
for ll = 1:length(pathDelays)
    hCoff = sqrt(1/2)*sqrt(averagePathPow(ll))*(randn + 1j*randn);    % 衰落信道
%     hCoff = 1;                                                        % 高斯信道
    p = raised_cosine(upSampRate,rollOffFactor,taud(ll),winLen);
%     genieChan(ll) = hCoff*max(p);
%     genieChan(ll) = hCoff;
    h = h + hCoff*(p.');
end
ht = h(winLen+1:end);
end

