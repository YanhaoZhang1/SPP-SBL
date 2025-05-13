% This demo shows how SPP-SBL algorithm reconstruct a block sparse chain signal
 

close all
clear all 
clc

iter = 1;
m = 130;
n = 512;
rand('seed',0);
randn('seed',0);

delta = 1e-1;
%%%%%%%%parameters of algorithm
p0hat = 0.8;  % p0hat=p(s_i=0)

num_iter_max = 5;

mu = 1e-7;

alpha = 0.99;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N=1;
gamma0=0;
gamma1=1;
p10 = 0.01;

sources=block_sparse(p0hat,p10,gamma0,gamma1,n);
%stem(sources);
SNR = 20;  


for it = 1:iter
 fprintf('\n\nRunning %d:\n',it);
    A=2*(rand(m,n)-.5); 
    % %%%%normalize columns
    A = A./(ones(m,1)*sqrt(sum(A.^2)));

    signal = A * sources;
    stdnoise = std(signal)*10^(-SNR/20);
    sign = stdnoise;
    noise = sign*randn(m,N);
    x=signal + noise;

    
   %% SPP-SBL

        c=10; d=1;
        beta = 1;
        tic;
        [Xspp]=SPP_SBL(x,A,stdnoise,beta,c,d); 
        t(it) = toc;
        mse_SPP(it)=(norm(sources - Xspp,'fro')/norm(sources,'fro'))^2;
        corr_SPP(it) = Xspp'*sources/(norm(Xspp)*norm(sources));

        fprintf('SPP_exact: time: %4.3f, MSE: %g\n',mean(t),mean(mse_SPP));
        fprintf('SPP_exact: time: %4.3f, Corr: %g\n',mean(t),mean(corr_SPP));

end


