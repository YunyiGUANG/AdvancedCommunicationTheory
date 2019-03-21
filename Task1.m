clear all
close all
clc

% read images
numofBits = 160*112*3*8;
[s1,x1,y1] = fImageSource('photo1.jpg',numofBits);
[s2,x2,y2] = fImageSource('photo2.jpg',numofBits);
[s3,x3,y3] = fImageSource('photo3.jpg',numofBits);

%% TASK 1.1 Transmitter
% modulation paramters
X = 7;
Y = 25;
coeffs1 = [1 0 0 1 1];
coeffs2 = [1 1 0 0 1];
phi = (X+2*Y)*pi/180; % phase shift
mseq1 = fMSeqGen(coeffs1);
mseq2 = fMSeqGen(coeffs2);
% find balanced gold sequence
d = 1 + mod((X+Y),12);
notbalanced=1;
while (notbalanced)
    goldseq1 = fGoldSeq(mseq1,mseq2,d);
    if(sum(goldseq1)~=-1)
        notbalanced=1;
    else
        notbalanced=0;
    end
    d = d+1;
end
goldseq2 = fGoldSeq(mseq1,mseq2,d);
goldseq3 = fGoldSeq(mseq1,mseq2,d+1);

s1=fDSQPSKModulator(s1,goldseq1,phi);
s2=fDSQPSKModulator(s2,goldseq2,phi);
s3=fDSQPSKModulator(s3,goldseq3,phi);

% channel parameters
paths = [1;1;1];
delay = [5;7;12];
beta = [0.4;0.7;0.2];
DOA = [30 0; 90 0; 150 0];
SNR = 40;
array = [0 0 0];
symbolsIn = [s1;s2;s3];
% transmit signals
symbolsOut = fChannel(paths,symbolsIn,delay,beta,DOA,SNR,array); 

%% TASK 1.2 Receiver
% Find delay
N = length(symbolsOut);
Nc = length(goldseq1); % length of the gold sequence
symbols = zeros(1,N-Nc);
for i=1:N-Nc
    seq = symbolsOut(i:i+Nc-1);
    symbols(i) = sum(seq.*goldseq1);
end
realpart = abs(symbols);
estimate_delay = mod(find(realpart==max(realpart),1),Nc)-1;

% Demodulation
symbolsIn = symbolsOut(estimate_delay+1:estimate_delay+length(s1));
demodulate_bits = fDSQPSKDemodulator(symbolsIn,phi,goldseq1);
fImageSink(demodulate_bits,numofBits,x1,y1);