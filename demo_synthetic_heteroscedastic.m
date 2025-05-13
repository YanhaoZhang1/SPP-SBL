% This demo shows how to use SPP-SBL to reconstruct the synthetic heteroscedastic signal.

close all
clear; clc;

rng(0,'twister')
M = 80;          % row number of the dictionary matrix 
N = 162;          % column number
blkNum = 4;       % nonzero block number
blkLen = 10;       % block length

SNR = 15;         % Signal-to-noise ratio
iterNum = 1;    % number of experiments


for it = 1: iterNum
    fprintf('\n\nRunning %d:\n',it);
    % Generate the known matrix with columns draw uniformly from the surface of a unit hypersphere
    Phi = randn(M,N);
    Phi = Phi./(ones(M,1)*sqrt(sum(Phi.^2)));

%%   Generate 'Heteroscedastic' data
    
    % Generate each block size randomly.
    minGap = 2; 
    partitionIndices = sort(randperm(blkLen*blkNum-1-(blkNum-1)*minGap, blkNum-1));
    partitionIndices = partitionIndices + (1:blkNum-1)*minGap;
    
    resultArray = diff([0 partitionIndices blkLen*blkNum]);
  
    if resultArray(end) <= minGap-1
        resultArray(end-1) = resultArray(end-1) - minGap+1;
        resultArray(end) = resultArray(end)+1;  
    end
    
    nums = resultArray;


    % Randomly select blkNum positions to place the blocks generate above.
    ind2 = zeros(1, blkNum);
    for i = 1:blkNum
        ind2(i) = randi([1, N-max(nums)+1]);
    end
    ind2 = sort(ind2);
    for i = 1:blkNum-1 
        sub = ind2(i+1)-ind2(i)+1;
        if sub<nums(i)
            ind2(i+1) = ind2(i+1)+(nums(i)-sub);
        end
    end


   % Generate blocks amplitude randomly
   rho=[];gen_B=[];gen_G=[];gen_Cov=[];blks={};
  
    for j=1:blkNum
        rho(j) = rand(1,1); 
        gen_B{j} = zeros(nums(j));
        eleB = 1; eleG = 1; 
        for i=1:(nums(j)-1)
            eleB = [eleB , rho(j)^i];
            eleG = [eleG , randi(20)];
        end
        gen_B{j} = toeplitz(eleB);
        gen_G{j} = diag(eleG);
        gen_Cov{j} = gen_G{j}*gen_B{j}*gen_G{j};

        % mean and covariance
        custom_mean = zeros(nums(j),1); 
        custom_covariance = gen_Cov{j}; 
        
        % multi-normal variable vector
        num_samples = 1; 
        random_numbers = mvnrnd(custom_mean, custom_covariance, num_samples);

        blks{j} = random_numbers;
    end
    
    % put blocks at locations
    Wgen = zeros(N,1);  blkLoc2 = [];
    for i = 1:blkNum
        nn = ind2(i)+nums(i)-1;
        mark = length(blkLoc2);
        blkLoc2 = [blkLoc2,ind2(i):(nn)];
        Wgen(ind2(i):(nn)) = blks{i};
    end
    if length(Wgen)>N
        Wgen=Wgen(1:end-(length(Wgen)-N));
    end

%% Noisy Measurements

    % noiseless signal
    signal = Phi * Wgen;
    
    % Observation noise   
    stdnoise = std(signal)*10^(-SNR/20);
    noise = randn(M,1) * stdnoise;

    % Noisy signal
    Y = signal + noise;
    IND = find(abs(Wgen)>0);
     
 
%% ================== SPP-SBL  ==============
    c=10; d=1;
    beta = 1;
    tic;
    [Xspp]=SPP_SBL(Y,Phi,stdnoise,beta,c,d); 
    t_exact_GS(it) = toc;
    mse_SPP(it)=(norm(Wgen - Xspp,'fro')/norm(Wgen,'fro'))^2;
    corr_SPP(it) = Xspp'*Wgen/(norm(Xspp)*norm(Wgen));
   
    fprintf('SPP_exact: time: %4.3f, MSE: %g\n',mean(t_exact_GS),mean(mse_SPP));
    fprintf('SPP_exact: time: %4.3f, Corr: %g\n',mean(t_exact_GS),mean(corr_SPP));


end

