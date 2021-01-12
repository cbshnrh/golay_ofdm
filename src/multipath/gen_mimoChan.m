function mimoChan = gen_mimoChan(pathDelays,averagePathGains,sampleRate,txNum,rxNum)
% pathDelays:多径的时延
% averagePathGains:多径的平均增益
% sampleRate:采样率 
winLen = 40;
mimoChan = zeros(winLen+1,rxNum,txNum);
for tt = 1:txNum
    for rr = 1:rxNum
        mimoChan(:,rr,tt) = gen_channel(pathDelays,averagePathGains,sampleRate,winLen);
    end
end
end

