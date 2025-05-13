 % This demo shows how to use SPP-SBL to reconstruct the synthetic multi-pattern signal.
close all
clear all 
clc
rng(0,'twister') 
iter = 1;

SNR = 10;    % Signal-to-noise ratio 
beta = 1;

n            = 162;           
m            = 80;
K            = 30;            % nonzero num
clustered_K  = 25;            % nonzero number in clusters
isolated_K   = K - clustered_K;% isolated num
blkNum       = 3;             % block num
cluster_lengths = [8, 4, 13];  % sum(cluster_lengths) = clustered_K
minGap       = 1;             

assert(sum(cluster_lengths)==clustered_K, 'Sum of cluster lengths must equal clustered_K');
assert(clustered_K + isolated_K == K,        'clustered_K + isolated_K must equal K');

signal = zeros(n,1);
used   = false(n,1);

ind2 = zeros(1, blkNum);
for j = 1:blkNum
    Bc = cluster_lengths(j);
    placed = false;
    while ~placed
        pos = randi([1, n - Bc + 1]);
        rngRange = pos:(pos+Bc-1);
        if all(~used(max(1,pos-minGap):min(n,pos+Bc-1+minGap)))
           
            rho = rand();
            eleB = rho.^(0:Bc-1);
            B    = toeplitz(eleB);
            G    = diag(randi(20,Bc,1));
            Cov  = G*B*G;
            blk  = mvnrnd(zeros(Bc,1), Cov, 1)'; 

            
            signal(rngRange) = blk;
            used(max(1,pos-minGap):min(n,pos+Bc-1+minGap)) = true;
            ind2(j) = pos;
            placed = true;
        end
    end
end


available = find(~used);
isolated_pos = [];
while numel(isolated_pos) < isolated_K
    idx = available(randi(numel(available)));
    
    if isempty(isolated_pos) || all(abs(isolated_pos - idx) > minGap)
        signal(idx) = 10*randn;  
        used(max(1,idx-minGap):min(n,idx+minGap)) = true;
        isolated_pos(end+1) = idx; 
        available = find(~used);
    end
end



for it = 1:iter
 fprintf('\n\nRunning %d:\n',it);
    % generate the measurement matrix
    Phi=randn(m,n);
    A=Phi./(ones(m,1)*sqrt(sum(Phi.^2)));
  
    stdnoise = std(signal)*10^(-SNR/20);
    sign = stdnoise;
    noise = sign*randn(m,1);
    x=A*signal + noise;

   %% SPP-SBL

    c=10; d=1;
    beta = 1;
    tic;
    [Xspp]=SPP_SBL(x,A,stdnoise,beta,c,d);
    t(it) = toc;
    mse_SPP(it)=(norm(signal - Xspp,'fro')/norm(signal,'fro'))^2;
    corr_SPP(it) = Xspp'*signal/(norm(Xspp)*norm(signal));
   
    fprintf('SPP_exact: time: %4.3f, MSE: %g\n',mean(t),mean(mse_SPP));
 

end







