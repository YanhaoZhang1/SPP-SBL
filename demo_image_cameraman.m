% Recovery of real world images using SPP-SBL.

clear;clc;close all;

rng(5,'twister') 
N=128;
% cameraman image (Set your own path here)
load('/Users/hhzhang/Desktop/SPP_SBL_public/cameraman.mat');
img=cameraman;

ss=[N,N]; % size of image
X0=imresize(img,ss);
X0=double(X0);
X=X0;
[a,b]=size(X);

% Discrete Wavelet Transform
load DWTM.mat  
M=64;
R=randn(M,a);
SNR=25; 
measure=R*X;
% Observation noise
stdnoise = std(measure)*10^(-SNR/20);
noise = randn(M,1) * stdnoise;
Y=measure+noise;
Phi=R*ww';

%% SPP-SBL

c=10; d=1;
beta = 1;
tic;
for i=1:N
    [Xspp]=SPP_SBL(Y(:,i),Phi,stdnoise(i),beta,c,d); % recover the im
    X_SPP(:,i) = Xspp;
end
t = toc;
X_SPP=ww'*X_SPP;  %  inverse-DWT transform
X_SPP=reshape(X_SPP,ss);
rnmse_SPP = sqrt(sum(sum((X_SPP-X0).^2,1),2)/sum(sum(X.^2,1),2));

fprintf('SPP_exact: time: %4.3f, MSE: %g\n',mean(t),mean(rnmse_SPP));

figure(1)
imshow(uint8(X_SPP));
title('SPP-SBL');


