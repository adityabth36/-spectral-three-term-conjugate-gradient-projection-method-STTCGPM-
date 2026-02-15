
% This is a demo of DFPM for solving constrained nonlinear weaker-monotone equations
% equations of the form
%   F(x)=0, x\in C, 
% where C is a nonempty closed convex set.

% model=1 means -F(zk)'*dk ≥ sigma*tk*norm(dk)^2
% model=2 means -F(zk)'*dk ≥ sigma*tk*norm(F(zk))*norm(dk)^2

clc;
clear all
close all

% set random number seed
rng(2025)

ITR_max = 2000;
% setup TXT document
fid_tex=fopen('mytext.txt','w'); 
% problem_set = 1:72;
 problem_set = [5:9 14:18 23:27 32:36 41:45 50:54 59:63 68:72 77:81];
% set parameters
np = length(problem_set); % from problem 1 to problem 8 % the number of the test problems
ns = 4;   % the number of the test algorithms
T = zeros(np,ns);
F = zeros(np,ns);
N = zeros(np,ns);
G = zeros(np,ns);

% set parameters for inertial derivative-free projection methods (IDFPMs)
% set parameters for STTCGPM
para1.Itr_max = ITR_max;
para1.gamma = 0.4;        % the initial guess
para1.sigma =  0.01;      % the coefficient of line search 
para1.tau = 0.3;         % the compression ratio
para1.alpha = 0.2;       % the coefficient of inertial step
para1.rho = 1.9;         % the relaxation factor 


% set parameters for inertial projection methods (FITTCGPM-PRP/FITTCGPM-DY/FITTCGPM-our)
para2.Itr_max = ITR_max;
para2.gamma = 0.4;         % the initial guess
para2.sigma = 0.01;         % the coefficient of line search 
para2.tau = 0.6;         % the compression ratio
para2.alpha = 0.1;       % the coefficient of inertial step
para2.rho = 1.8;         % the relaxation factor 

% set parameters for IMSMNE
para3.Itr_max = ITR_max;
para3.gamma = 1;         % the initial guess
para3.sigma = 0.001;         % the coefficient of line search 
para3.tau = 0.36;         % the compression ratio
para3.alpha = 0.35;       % the coefficient of inertial step
para3.rho = 1.87;         % the relaxation factor 

% run
for index=1:40  %np
    Num = problem_set(index);
    [name,n] = init(Num);
    progress_r = [];
    for repeats = 1:5
        %  if (Num<=9 && Num>=5)
        %     x0 = 2*rand(n,1)-1;
        % else
            x0 = 2*rand(n,1);    % x0\in (0,10). x0 = 4*rand(n,1)
        % end
        [T1,NFF1,NI1,G1] = DFPM(Num,'STTCGPM',1,2,x0,para1); % acceleration
        [T2,NFF2,NI2,G2] = DFPM(Num,'FITTCGPM-PRP',2,1,x0,para2);
        [T3,NFF3,NI3,G3] = DFPM(Num,'FITTCGPM-DY',2,1,x0,para2);
        [T4,NFF4,NI4,G4] = DFPM(Num,'IMSMNE',3,2,x0,para3);
        progress_r = [progress_r;NI1,NFF1,T1,G1,NI2,NFF2,T2,G2,NI3,NFF3,T3,G3,NI4,NFF4,T4,G4];%,NI5,NFF5,T5,G5];
    end
    TM = mean(progress_r);
    fprintf(fid_tex,'%s %d & %.1f/%.1f/%.3f/%.2e & %.1f/%.1f/%.3f/%.2e\n& %.1f/%.1f/%.3f/%.2e & %.1f/%.1f/%.3f/%.2e\\\\ \r\n', ... 
                name,n,TM);
    T(index,:) = [TM(3),TM(7),TM(11),TM(15)];
    F(index,:) = [TM(2),TM(6),TM(10),TM(14)];
    N(index,:) = [TM(1),TM(5),TM(9),TM(13)];
    G(index,:) = [TM(4),TM(8),TM(12),TM(16)];
end
%% close file
fclose(fid_tex);

%% Draw a picture
clf;   %Clear Figure
    
figure(1);
%subplot(2,2,1);
perf(T,'logplot');
%(Time performance)
%set(gca,'ylim',[0.3,1]);
xlabel('\tau','Interpreter','tex');
ylabel('\rho(\tau)','Interpreter','tex');
legend('STTCGPM','FITTCGPM-PRP','FITTCGPM-DY','IMSMNE','Location','southeast');
% %subplot(2,2,2);
figure(2);
perf(F,'logplot');
%('Objective function calculation performance');
% set(gca,'ylim',[0.1,1]);
xlabel('\tau','Interpreter','tex');                    
ylabel('\rho(\tau)','Interpreter','tex');               
legend('STTCGPM','FITTCGPM-PRP','FITTCGPM-DY','IMSMNE','Location','southeast');
%subplot(2,2,3);
figure(3);
perf(N,'logplot');
%('Iteration number performance');
%set(gca,'ylim',[0.5,1]);
xlabel('\tau','Interpreter','tex');
ylabel('\rho(\tau)','Interpreter','tex');
legend('STTCGPM','FITTCGPM-PRP','FITTCGPM-DY','IMSMNE','Location','southeast');
