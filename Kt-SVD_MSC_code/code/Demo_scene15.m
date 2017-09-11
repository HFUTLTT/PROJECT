clear;
addpath('tSVD','proxFunctions','solvers','twist');
addpath('ClusteringMeasure', 'LRR', 'Nuclear_norm_l21_Algorithm', 'unlocbox');

load('data/scene15.mat');                          
cls_num = length(unique(gt));
%% Note: each column is an sample (same as in LRR)
%% 
%data preparation...
 X{1} = X1; X{2} = X2; X{3} = X3;
 K = length(X); N = size(X{1},2); %sample number
 for v=1:3
     [X{v}]=NormalizeData(X{v}); 
 end
Kx{3}=LinearKernel(X{3});
Kx{2}=GaussKernel(X{2});
Kx{1}=LinearKernel(X{1});
 for v=1:3
     [V,S,Vt]=svd(Kx{v},'econ'); 
     Rk{v}=length(diag(S));
     Sk{v}=sqrt(diag(S));
     Sk{v}=Sk{v}(1:Rk{v});
     Vk{v}=V(:,1:Rk{v});
 end
% Initialize...

for k=1:K
    Z{k} = zeros(N,N); %Z{2} = zeros(N,N);
    W{k} = zeros(N,N);
    G{k} = zeros(N,N);
    %E{k} = zeros(size(X{k},1),N); %E{2} = zeros(size(X{k},1),N);
    P{k} = eye(N,N);
    %Y{k} = zeros(size(X{k},1),N); %Y{2} = zeros(size(X{k},1),N);
    Y{k} = zeros(N,N);
end

w = zeros(N*N*K,1);
g = zeros(N*N*K,1);
dim1 = N;dim2 = N;dim3 = K;
myNorm = 'tSVD_1';
sX = [N, N, K];
%set Default
%% 
parOP         =    false;
ABSTOL        =    1e-6;
RELTOL        =    1e-4;
%%
Isconverg = 0;epson = 1e-7;
lambda =0.0018; %1.5 best
iter = 0;
mu = 10e-5; max_mu = 10e13; pho_mu =2;
rho = 10e-5; max_rho = 10e13; pho_rho = 2;
tic;

while(Isconverg == 0)
    fprintf('----processing iter %d--------\n', iter+1);
    for k=1:K
        %1 update Z^k
        Z{k}=(mu*eye(N,N)+Y{k}+rho*G{k}-mu*P{k}-W{k})/(rho+mu);
       %%        
        %2 update P^k
               
        C{k}=eye(N,N)-Z{k}+Y{k}/mu;
        tau=mu/lambda;
        tu=Vk{k}'*C{k};
        % upadte Pki
        for i=1:N
            tui=tu(:,i);
            ci=C{k}(:,i);
            condition=norm(tui./Sk{k},2);
            if condition>1/tau
             %% 
               %alpha=get_alpha_m(tui,Sk{k},tau);
               alpha=get_alpha(tui,Sk{k},tau,Rk{k});
             %% 
               tmp1=zeros(Rk{k},1);
               for j=1:Rk{k}
                  tmp1(j)=Sk{k}(j)^2*tui(j)/(alpha+Sk{k}(j)^2);
               end
 
               P{k}(:,i)=ci-Vk{k}*tmp1; 
            else
               P{k}(:,i)= ci-Vk{k}*tui; 
            end
         end
        %3 update Yk
        Y{k} = Y{k} + mu*(eye(N,N)-Z{k}-P{k});
    end   
    %4 update G
    Z_tensor = cat(3, Z{:,:});
    W_tensor = cat(3, W{:,:});
    z = Z_tensor(:);
    w = W_tensor(:);    
    %twist-version
   [g, objV] = wshrinkObj(z + 1/rho*w,1/rho,sX,0,3)   ;
    G_tensor = reshape(g, sX);    
    %5 update W
    w = w + rho*(z - g);
    
    %record the iteration information
    history.objval(iter+1)   =  objV;

    %% coverge condition
    Isconverg = 1;
    for k=1:K
        if (norm(eye(N,N)-Z{k}-P{k},inf)>epson)
            history.norm_Z = norm(eye(N,N)-Z{k}-P{k},inf);
            fprintf('    norm_Z %7.10f    ', history.norm_Z);
            Isconverg = 0;
        end
        
        G{k} = G_tensor(:,:,k);
        W_tensor = reshape(w, sX);
        W{k} = W_tensor(:,:,k);
        if (norm(Z{k}-G{k},inf)>epson)
            history.norm_Z_G = norm(Z{k}-G{k},inf);
            fprintf('norm_Z_G %7.10f    \n', history.norm_Z_G);
            Isconverg = 0;
        end
    end
   
    if (iter>200)
        Isconverg  = 1;
    end
    iter = iter + 1;
    mu = min(mu*pho_mu, max_mu);
    rho = min(rho*pho_rho, max_rho);
end
S = 0;
for k=1:K
    S = S + abs(Z{k})+abs(Z{k}');
end
C = SpectralClustering(S,cls_num);
[A nmi avgent] = compute_nmi(gt,C);
ACC = Accuracy(C,double(gt));
[f,p,r] = compute_f(gt,C);
[AR,RI,MI,HI]=RandIndex(gt,C);
toc;
%save('my_new_COIL20MV_res.mat','S','ACC','nmi','AR','f','p','r');




