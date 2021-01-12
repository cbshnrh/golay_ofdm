function ht = gen_channel(pathDelays,averagePathGains,sampleRate,winLen)
% pathDelays:�ྶ��ʱ��
% averagePathGains:�ྶ��ƽ������
% sampleRate:������  �˴�sampleRate=1(���ù�)
% winLen:sinc�����Ĵ��ڴ�С
% ht:���ɵ�ʱ��ྶ�ŵ�(����������)
% winLen = 40;
averagePathPow = db2pow(averagePathGains);
taud = pathDelays*sampleRate;
rollOffFactor = 0.8;
upSampRate = 1;
h = zeros(1,winLen*2*upSampRate+1);
% genieChan = zeros(1,length(pathDelays));
for ll = 1:length(pathDelays)
    hCoff = sqrt(1/2)*sqrt(averagePathPow(ll))*(randn + 1j*randn);    % ˥���ŵ�
%     hCoff = 1;                                                        % ��˹�ŵ�
    p = raised_cosine(upSampRate,rollOffFactor,taud(ll),winLen);
%     genieChan(ll) = hCoff*max(p);
%     genieChan(ll) = hCoff;
    h = h + hCoff*(p.');
end
ht = h(winLen+1:end);
end

