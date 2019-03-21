clear all
close all
clc

% read signals
numofBits = 160*112*3*8;
[s1,x1,y1] = fImageSource('photo1.jpg',numofBits);
[s2,x2,y2] = fImageSource('photo2.jpg',numofBits);
[s3,x3,y3] = fImageSource('photo3.jpg',numofBits);

%% Transmitter
% Modulation parameters
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

s1=fDSQPSKModulator(s1,goldseq1,phi);
s2=fDSQPSKModulator(s2,goldseq2,phi);
s3=fDSQPSKModulator(s3,goldseq3,phi);

% uniform circular array with 5 isotropic antennas
array = zeros(5,3);
l = 1/sqrt(2*(1-cos(72*pi/180)));
for i=1:5
    angles = 30 + 360*(i-1)/5;
    pha = angles*pi/180;
    array(i,:) = l*[cos(pha) sin(pha) 0];
end
% channel parameters
paths = [1;1;1];
delay = [5;7;12];
beta = [0.4;0.7;0.8];
DOA = [30 0; 90 0; 150 0];
SNR = 40;
symbolsIn = [s1;s2;s3];
% transmit signals
symbolsOut = fChannel(paths,symbolsIn,delay,beta,DOA,SNR,array);

% channel estimation
user = 1; % desired user
numofPaths = paths(user); % number of paths of the desired user
[delay_estimate, DOA_estimate, beta_estimate]=fChannelEstimation(symbolsOut,goldseq1,numofPaths);

%% Receiver
% WIENER-HOPF Beamformer
Sd = spv(array,DOA_estimate);
Rxx = symbolsOut*symbolsOut';
alpha = 1;
wopt = alpha*inv(Rxx)*Sd;
symbols = wopt'*symbolsOut;

% Demodulation
symbolsIn = symbols(delay_estimate+1:delay_estimate+length(s1));
demodulate_bits = fDSQPSKDemodulator(symbolsIn,phi,goldseq1);
fImageSink(demodulate_bits,numofBits,x1,y1);