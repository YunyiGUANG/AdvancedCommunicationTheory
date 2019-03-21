clear all
close all
clc

% read images
numofBits = 160*112*3*8;
[s1,x1,y1] = fImageSource('photo1.jpg',numofBits);
[s2,x2,y2] = fImageSource('photo2.jpg',numofBits);
[s3,x3,y3] = fImageSource('photo3.jpg',numofBits);

%% Transmitter
% modulation parameters
X = 7;
Y = 25;
coeffs1 = [1 0 0 1 1];
coeffs2 = [1 1 0 0 1];
phi = (X+2*Y)*pi/180;
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

% modulation
s1=fDSQPSKModulator(s1,goldseq1,phi);
s2=fDSQPSKModulator(s2,goldseq2,phi);
s3=fDSQPSKModulator(s3,goldseq3,phi);

% channel parameters
paths = [3;1;1];
delay = [mod(X+Y,4);4+mod(X+Y,5);9+mod(X+Y,6);8;13];
beta = [0.8;0.4*exp(-1j*40*pi/180);0.8*exp(1j*80*pi/180);0.5;0.2];
DOA = [30 0; 45 0; 20 0; 80 0; 150 0];
SNR = 40;
array = [0 0 0];
symbolsIn = [s1;s2;s3];
% transmit signals
symbolsOut = fChannel(paths,symbolsIn,delay,beta,DOA,SNR,array);

%% Receiver
% find delay
user = 1; % the desired user
N = length(symbolsOut);
Nc = length(goldseq1); % the length of gold sequence
symbols = zeros(1,N-Nc);
for i=1:N-Nc
    seq = symbolsOut(i:i+Nc-1);
    symbols(i) = sum(seq.*goldseq1);
end
realpart = abs(symbols);
realpart = reshape(realpart,[15,length(realpart)/15]);
realpart = sum(realpart.');
[~,idx] = sort(realpart,'descend');
numofPaths = paths(user); % the number of paths of the desired user
estimate_delay = idx(1:numofPaths)-1; % find the index of the largest N value( N=numofPaths)
estimate_delay = sort(estimate_delay);

% Demodulation
symbolsIn = zeros(1,length(s1));
% combine all the paths of desired user
for i=1:length(estimate_delay)
    symbolsIn = symbolsIn + symbolsOut(estimate_delay(i)+1:estimate_delay(i)+length(s1))/beta(i);
end
demodulate_bits = fDSQPSKDemodulator(symbolsIn,phi,goldseq1);
fImageSink(demodulate_bits,numofBits,x1,y1);