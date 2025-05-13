% This demo shows how to use SPP-SBL to reconstruct real-world audio signal.

close all
clear; clc;
rng(0,'twister') 

info =audioinfo('-0SdAVK79lg.wav');
[audio,Fs] = audioread('-0SdAVK79lg.wav');


audiolength = 480;
t = 1:1:audiolength;
audio = audio(:,1);
audio = audio(t);
figure(1),
plot(t,audio);
xlabel('Time');
ylabel('Audio Signal');

rng(0,'twister') 
M = 150; % measurement number
N = audiolength;
SNR = 15;  

Wgen = dct(audio);
iterNum = 100;    % number of experiments


for it = 1: iterNum

    fprintf('\n\nRunning %d:\n',it);
    
    % Generate the known matrix with columns draw uniformly from the surface of a unit hypersphere
    Phi = randn(M,N);
    Phi = Phi./(ones(M,1)*sqrt(sum(Phi.^2)));
    
    
    %% Measurements
    % noiseless signal
    signal = Phi * Wgen;
    
    % Observation noise   
    stdnoise = std(signal)*10^(-SNR/20);
    noise = randn(M,1) * stdnoise;
    
    % Noisy signal
    Y = signal + noise; 
    IND = find(abs(Wgen)>0);

    %% SPP-SBL
    c=10; d=1;
    beta = 1; % initial
    tic;
    [Xspp]=SPP_SBL(Y,Phi,stdnoise,beta,c,d); % recover the im
    t(it) = toc;
    mse_SPP(it)=(norm(Wgen - Xspp,'fro')/norm(Wgen,'fro'))^2;
    fprintf('SPP_exact: time: %4.3f, MSE: %g\n',mean(t),mean(mse_SPP));

end




