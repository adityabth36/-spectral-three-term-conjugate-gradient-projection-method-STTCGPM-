
% consider the regularized decentralized logistic regression
%
% $$ \displaystyle\min_x \frac{1}{m}\sum_{i=1}^m \log(1+ \exp(-b_ia_i^Tx)) 
% + \lambda \|x\|^2/2,$$
% 
% where $(a_i,b_i)_{i=1}^m$




% model=1 means -F(zk)'*dk ?? sigma*tk*norm(dk)^2
% model=2 means -F(zk)'*dk ?? sigma*tk*norm(F(zk))*norm(dk)^2
% model=3 means -F(zk)'*dk ?? sigma*tk*norm(F(zk))/(1+norm(F(zk)))*norm(dk)^2
% model=4 means -F(zk)'*dk ?? sigma*tk*max(lambda,min(nu,norm(Fz_new,2)))*norm(dk)^2

clc;
clear all
close all

%%
% set random seed 
clear;
seed = 97006855;
ss = RandStream('mt19937ar','Seed',seed);
RandStream.setGlobalStream(ss);

ITR_max = 2000;

% set parameters for HCGPM
para1.Itr_max = ITR_max;
para1.gamma = 0.4;         % the initial guess
para1.sigma = 0.01;      % the coefficient of line search 
para1.tau = 0.3;         % the compression ratio
para1.rho = 1.9;         % the relaxation factor 
% para2.mu =10; 
% para2.lambda_k=15;
para1.alpha = 0.2;       % the coefficient of inertial step
% para2.eps = 10^-6;

% set parameters for MITTCGP
para2.Itr_max = ITR_max;
para2.gamma = 0.4;         % the initial guess
para2.sigma = 0.01;      % the coefficient of line search 
para2.tau = 0.6;         % the compression ratio
para2.rho = 1.8;         % the relaxation factorsigma=0.01;  % sigma=0.01
para2.alpha = 0.1;       % the coefficient of inertial step

% set parameters for IMSMNE
para3.Itr_max = ITR_max;
para3.gamma = 1;         % the initial guess
para3.sigma = 0.001;         % the coefficient of line search 
para3.tau = 0.36;         % the compression ratio
para3.alpha = 0.35;       % the coefficient of inertial step
para3.rho = 1.87;         % the relaxation factor 

lambda = 1e-1;
% setup TXT document
fid_tex=fopen('mytext.txt','w'); 

for i=10:11
    dataset = {'a1a.t','a2a.t','a3a.t','a4a.t','a5a.t','a6a.t','a7a.t','a8a.t','a9a.t','ionosphere_scale.t',...
    'w1a.t','w2a.t','w3a.t','w4a.t','w5a.t','w6a.t','w7a.t','w8a.t'};
    [b,A] = libsvmread(dataset{i});
    [m,n] = size(A);
    fprintf('name=%s, m=%d, n=%d, lambda=%.2f\n',dataset{i},m,n,lambda);

    % create function handle
    fun = @(x) lr_loss(A,b,m,x,lambda);
    progress_r = [];
    for repeat=1:5
        % set the initial point
        x0 = 4*(rand(n,1)-0.5); %4*(rand(n,1)-0.5);    % 2*(rand(n,1)-0.5)
        x_1 = x0;   %4*(rand(n,1)-0.5); 
       %% start comparison %%
        [T1,NFF1,NI1,G1] = IDFPM(fun,'STTCGPM',1,2,x0,x_1,para1);
        [T2,NFF2,NI2,G2] = IDFPM(fun,'FITTCGPM-PRP',2,1,x0,x_1,para2);
        [T3,NFF3,NI3,G3] = IDFPM(fun,'FITTCGPM-DY',2,1,x0,x_1,para2);
        [T4,NFF4,NI4,G4] = IDFPM(fun,'IMSMNE',1,2,x0,x_1,para3);
        %InerDFPI(fun,'IDFPI',3,1,x0,x_1,para8); 
%         ISTCP(fun,'ISTCP',2,2,x0,x_1,para5);
        progress_r = [progress_r;NI1,NFF1,T1,G1,NI2,NFF2,T2,G2,NI3,NFF3,T3,G3,NI4,NFF4,T4,G4];
    end
    TM = mean(progress_r);
    fprintf(fid_tex,'%s & %d & %d & %.1f/%.1f/%.3f/%.2e & %.1f/%.1f/%.3f/%.2e & %.1f/%.1f/%.3f/%.2e & %.1f/%.1f/%.3f/%.2e\\\\\\ \r\n',...
                dataset{i},m,n,TM);%NI1,NFF1,T1,G1,NI2,NFF2,T2,G2); 
%     fprintf(fid_tex,'%s & %d & %d & %.1f/%.1f/%.3f/%.2e & %.1f/%.1f/%.3f/%.2e\n& %.1f/%.1f/%.3f/%.2e & %.1f/%.1f/%.3f/%.2e\\\\ \r\n',...
%                 dataset{i},m,n,NI1,NFF1,T1,G1,NI2,NFF2,T2,G2,NI3,NFF3,T3,G3,NI4,NFF4,T4,G4); 

end
%% close file 
fclose(fid_tex);