clear all
close all
clc

% load files
load('yg1618.mat');

% known parameters
coeff1 = [1 0 0 1 0 1];
coeff2 = [1 0 1 1 1 1];
mseq1 = fMSeqGen(coeff1);
mseq2 = fMSeqGen(coeff2);
goldseq = fGoldSeq(mseq1,mseq2,phase_shift);

% channel estimation
numofPaths = 3;
[delay_estimate, DOA_estimate, beta_estimate]=fChannelEstimation(Xmatrix,goldseq,numofPaths);

Nc = length(goldseq);
array = zeros(5,3);
% uniform circular array with 5 isotropic antennas
l = 1/sqrt(2*(1-cos(72*pi/180)));
for i=1:5
    angles = 30 + 360*(i-1)/5;
    pha = angles*pi/180;
    array(i,:) = l*[cos(pha) sin(pha) 0];
end

% STAR Algorithm
C = [goldseq';zeros(Nc,1)];
J = [zeros(2*Nc-1,1)' 0; eye(2*Nc-1) zeros(2*Nc-1,1)]; % shifting matrix
for i=1:3
    Jc = J^delay_estimate(i)*C;
    Sd = spv(array,DOA_estimate(i,:));
    h(i,:) = kron(Sd,Jc); % STAR manifold vector
end
H = [h(1,:).',h(2,:).',h(3,:).']; % STAR manifold matrix
[M,N] = size(Xmatrix);
X = zeros(2*Nc*M, N/Nc-1);
% transposition and vectorisation of data Xmatrix
for i=1:M
    for j=1:N/Nc-1
        T = Xmatrix(i,(j-1)*Nc+1:(j+1)*Nc);
        X((i-1)*2*Nc+1:i*2*Nc,j) = T.';
    end
end
% spatiotemoral-RAKE
w = H*Beta_1;
symbols = w'*X;

% demodulation
bitsOut = fDSQPSKDemodulator(symbols,phi_mod*pi/180);
demodulate_bits = reshape(bitsOut,[8,length(bitsOut)/8]);
bitsOut = bi2de(demodulate_bits','left-msb');
charOut = char(bitsOut.')