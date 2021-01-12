% test mimo chan
txNum = 4;
rxNum = 1;
N = 1000;
txSig = 1/sqrt(2)*(randn(txNum,N));

pathDelays = [0 0.042 0.101 0.129 0.149 0.245 0.312 0.410 0.469 0.528]*1e-6;  
averagePathGains = [-5.2 -6.4 -8.4 -9.3 -10 -13.1 -15.3 -18.5 -20.4 -22.4];   
sampleRate = 61.44*1e6;   

mimoChan = gen_mimoChan(pathDelays,averagePathGains,sampleRate,txNum,rxNum); 
rxSig = mimo_channel(mimoChan,txSig.'); 