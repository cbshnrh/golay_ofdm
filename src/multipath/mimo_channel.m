function rxSig = mimo_channel(mimoChan,txSig)
% mimoChan��MIMO�ŵ�ʱ����Ӧ 
% txSig�������ź�
[~,rxNum,txNum] = size(mimoChan);
[sigLen,~] = size(txSig);
rxSig = zeros(sigLen,rxNum);
for rr = 1:rxNum
    for tt = 1:txNum
        ht = mimoChan(:,rr,tt);
        y = conv(ht,txSig(:,tt));
        rxSig(:,rr) = rxSig(:,rr) + y(1:sigLen);
    end
end
end

