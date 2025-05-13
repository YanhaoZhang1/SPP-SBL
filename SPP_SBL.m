function [x_new,eta]=SPP_SBL(y,A,sigma,beta,c,d)

%%
% Input: 
% y: measurements;
% A: sensing matrix;
% sigma: square root of the noise covariance;
% beta: parameter controling the relevance between the elements;
% Output:
% x_new: estimated sparse signal
% eta: estimated coupling parameters in each iteration


[m,n]=size(A);
eta = [];
iter=0;
iter_mx=100;
beta = beta*ones(n,1);
beta(n) = 0;
betal = [beta(2:n);0];
betar = [0;beta(1:n-1)];

D=eye(n);
sigma2=sigma^2;
alpha_new=ones(n,1);
var_new=inv(A'*A/sigma2+D);
mu_old=ones(n,1);
mu_new=1/sigma2*var_new*A'*y;
a=0.5;
b=1e-10;

root = zeros(n,1);


while iter<iter_mx& norm(mu_new-mu_old)>1e-6

    iter=iter+1;
    mu_old=mu_new;
    mul=[mu_new(2:n);0];
    mur=[0;mu_new(1:n-1)];
    var=diag(var_new);
    varl=[var(2:n);0];
    varr=[0;var(1:n-1)];

    % update alpha,beta
    for i = 1:n
        E=mu_new(i)^2+beta(i)*mul(i)^2+betar(i)*mur(i)^2+var(i)+beta(i)*varl(i)+betar(i)*varr(i);
        alpha_new(i) = a/(0.5*E+b);
         % force-zero
        idx1=find(alpha_new>1e10);
        alpha_new(idx1)=1e10;
        alf = [alpha_new(2:n); 0];                                %   left-shifted version
        arf=[0; alpha_new(1:n-1)];                                  %   right-shifted version
        allf = [alpha_new(3:n); 0;0];

        if i==n
            root(i)=0;
            beta(i) = root(i);
            betal = [beta(2:n);0];
            betar = [0;beta(1:n-1)];
        else
            % Div-power prior for eta
            B = 0.5*(alf(i)*(mu_new(i)^2+var(i))+alpha_new(i)*(mul(i)^2+varl(i)));
            C = alpha_new(i)+betar(i)*arf(i);
            EE = alf(i)+betal(i)*allf(i);
        
            M = 2*(B+d)*alpha_new(i)*alf(i);
            N = 2*(B+d)*(C*alpha_new(i)+alf(i)*EE)-2*c*alpha_new(i)*alf(i);
            Q = 2*(B+d)*C*EE+(1-2*c)*(C*alpha_new(i)+alf(i)*EE);
            R = 2*(1-c)*C*EE;
    
            %function
            p = [];
            p = [M,N,Q,R];
            
            solution = roots(p);
              
            real_roots = solution(imag(solution) == 0);

            if isempty(real_roots)
                fprintf('No real solution');
            else
                root(i) = max(real_roots);
        %         disp(['The max real solution is: ', num2str(max_real)]);
            end
            solution(solution==root(i))=[];
            % record_root{iter,i,:} = solution;
            beta(i) = root(i);
            betal = [beta(2:n);0];
            betar = [0;beta(1:n-1)];
        end
        
    end

    D=diag(alpha_new+betar.*arf+beta.*alf);
    var_new=inv(A'*A/sigma2+D);
    mu_new=1/sigma2*var_new*A'*y;
    eta(:,iter) = beta;

end
x_new=mu_new;
