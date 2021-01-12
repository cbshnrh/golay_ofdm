numMonte = 1e3;
snrDbSet = 8:1:20;
sig;
for nn = 1:numMonte
    h = gen_chan();
    ns = (randn(1,sigLen)+1j*randn(1,sigLen))/sqrt(2);
    for snrIndx = 1:length(snrDbSet)
        if numErrBit(snrIndx)>=100
            continue
        end
        snrDb = snrDbSet(snrIndx);
        x = h*sig+ns*db2mag(snrDb);
    end
        
    end2 