function mimoChan = gen_mimoChan(pathDelays,averagePathGains,sampleRate,txNum,rxNum)
% pathDelays:�ྶ��ʱ��
% averagePathGains:�ྶ��ƽ������
% sampleRate:������ 
winLen = 40;
mimoChan = zeros(winLen+1,rxNum,txNum);
for tt = 1:txNum
    for rr = 1:rxNum
        mimoChan(:,rr,tt) = gen_channel(pathDelays,averagePathGains,sampleRate,winLen);
    end
end
end

