%the inertial derivative-free projection method (IDFPM) for solving 
%% constrained nonlinear weaker-monotone equations of the form
%   F(x)=0, x\in K, 
% where K is a nonempty closed convex set.

function [Tcpu,NF,Itr,NormF] = DFPM(NO,method,Switch,model,x0,para) 
format long
% start the clock
tic;

%% Initial information
[nprob,~]=init(NO);
%% the stopping criterion
epsilon=1e-6;
epsilon1=1e-7;

%% the line search parameters and relaxation factor
k_max = para.Itr_max;   % the maximum number of iterations
gamma = para.gamma;     % the initial guess
sigma = para.sigma;     % the coefficient of line search 
tau = para.tau;         % the compression ratio
alpha = para.alpha;     % the coefficient of the inertial step
rho = para.rho;         % the relaxation factor 
% alpha_try = 0.1;

% eps = para.eps;         % Truncation creterion
 fprintf('%s & %s & LSmodel=%d & gamma=%.4f & sigma=%.4f & tau=%.4f & Switch=%.4f & rho=%.4f\n', ... 
     nprob,method,model,gamma,sigma,tau,Switch,rho);
%% compute the search direction
Fx0 = feval(nprob,x0,1);   % evaluate the function value specified by nprob at x0
NF = 1;  
NormFx0 = norm(Fx0);                   
x0_old = x0;
L1 = 0;
     
for k=1:k_max
    
    if k==1 && NormFx0<=epsilon
        L1 = 1;
        NormF = NormFx0; % the final norm of equations
        break; 
    end
    ell = norm(x0-x0_old);

    %% 
    p=1.9; q1=0.2; q2=1;
     if k==1
        alpha_k = alpha;
    else
        if Switch==1
            if  ell <=1
            alpha_k = min(alpha,1/(k^p * ell^(q1)));
            else 
            alpha_k = min(alpha,1/(k^p * ell^(q2)));
            end
        elseif Switch==2
            alpha_k = min(alpha,1/(k*ell)^2);
        elseif Switch==3
           alpha_k = min(alpha,1/(k^2 * ell^(1)));
        else
            alpha_k = alpha;
        end
    end
    %% compute the inertial step %%
    v0 = x0+alpha_k*(x0-x0_old);
    Fv0 = feval(nprob,v0,1);
    NF = NF+1;
    NormFv0 = norm(Fv0);
    if NormFv0<=epsilon
        L1 = 1;
        NormF = NormFv0;   % the final norm of equations
        break; 
    end
    
    %% compute the initial direction %%
    if k==1
        dk = -Fv0;
    else
        % update the search direction
        switch method   
           case 'STTCGPM'
               mu =10;
               lambda_k=15;
               w0 = Fv0-Fy0_old;
               NormFv02 = NormFv0^2 ;
               U_k = mu *(Normdk2 + NormFv02+ abs(dk'*w0));
               Fv0tdk = Fv0'*dk;
               beta_k = NormFv02/U_k - ((NormFv02 * (Fv0tdk))/ U_k^2);
               theta_k= lambda_k+beta_k*((Fv0tdk)/(NormFv02));
               p_k = Fv0 ;
               Normpk =NormFv0;
               dk=-theta_k.*Fv0 + beta_k .*dk + ((Fv0tdk)/(Normdk*Normpk))*p_k;
           case 'FITTCGPM-PRP'
                 w0 = Fv0-Fy0_old;
                 NormFy0 = norm(Fy0_old);
                 NormFy02 = norm(Fy0_old)^2;
                 betak=Fv0'*w0/( NormFy02);
                 nuk=NormFv0/max(NormFy0,(abs(betak))*Normdk);
                 thetak=-0.055*Fv0'*Fy0_old/ NormFy02;   
                 dk=-Fv0+0.001*nuk*betak*dk+thetak*Fy0_old;
           case'FITTCGPM-DY'
                 w0 = Fv0-Fy0_old;
                 NormFy0 = norm(Fy0_old);
                 betak=NormFv0^2/(dk'*w0);
                 nuk=NormFv0/max(NormFy0,(abs(betak))*Normdk);
                 thetak=-0.055*Fv0'*Fy0_old/NormFy0^2;   
                 dk=-Fv0+0.001*nuk*betak*dk+thetak*Fy0_old;   
          case 'FITTCGPM-NEW'
                w0 = Fv0-Fy0_old;
                U_k = mu *(Normdk2 + NormFv0^2+ abs(dk'*w0));
                betak = NormFv0^2/U_k - ((NormFv0^2 * (Fv0'*dk))/ U_k^2);
                 nuk=norm(Fv0)/max(norm(Fy0_old),(abs(betak))*Normdk);
                thetak=-0.055*Fv0'*Fy0_old/norm(Fy0_old)^2;  
                dk=-Fv0+0.001*nuk*betak*dk+thetak*Fy0_old;
           case 'IMSMNE'
                t1=0.3; t2=10^-3; D=10; r0=3; tauk=6.2 ; 
                w0 = Fv0-Fy0_old;
                sk=v0-y0_old;
                Normsk =sk'*sk;
                h_k = w0 + D*norm(Fy0_old)^r0 * sk + max(0,((-sk'*w0)/(Normsk)))*sk;
                Normhk =norm(h_k);
                k_a = max(tauk*Normdk*Normhk, dk'*h_k);
                betak= (Fv0'*h_k)/(k_a) - ((Normhk^2 * Fv0'*dk)/(k_a^2))  ;
                theta = ((sk'*h_k)/(Normsk));
                theta_k =  min(max(theta, t1), t2);
                dk = - theta_k * Fv0 + betak * dk ;     
            otherwise
                disp('Input error! Please check the input method');
        end
    end
    Normdk = norm(dk);
    if Normdk<epsilon1
        L1 = 1;
        NormF = NormFv0;
        break;
    end
    Normdk2 = Normdk^2;
    Fy0_old = Fv0;
    
    %%% Start Armijo-type line search  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % model=1 means -F(zk)'*dk ≥ sigma*tk*norm(dk)^2
    % model=2 means -F(zk)'*dk ≥ sigma*tk*norm(F(zk))*norm(dk)^2
    
    if model==1
        t = gamma;
        z0_new = v0+t*dk;
        Fz0_new = feval(nprob,z0_new,1);
        NF = NF+1;
        Fz0_newtdk = -Fz0_new'*dk;
        % check the Armijo-type line search condition
        while Fz0_newtdk < sigma*t*Normdk2 %&& t>10^-6  
            % the Armijo-type line search condition violated
            t = t*tau;
            z0_new = v0+t*dk;
            Fz0_new = feval(nprob,z0_new,1);
            NF = NF+1;
            Fz0_newtdk = -Fz0_new'*dk;
        end %%% End Armijo-type line search %%%
        NormFz0_new = norm(Fz0_new);
    else
        t = gamma;
        z0_new = v0+t*dk;
        Fz0_new = feval(nprob,z0_new,1);
        NF = NF+1;
        NormFz0_new = norm(Fz0_new);
        Fz0_newtdk = -Fz0_new'*dk;
        % check the Armijo-type line search condition
        while Fz0_newtdk < sigma*t*NormFz0_new*Normdk2 && t>10^-6  
            % the Armijo-type line search condition violated
            t = t*tau;
            z0_new = v0+t*dk;
            Fz0_new = feval(nprob,z0_new,1);
            NF = NF+1;
            NormFz0_new = norm(Fz0_new);
            Fz0_newtdk = -Fz0_new'*dk;
        end %%% End Armijo-type line search %%%
    end 
    Fz0 = Fz0_new;
    NormFz0 = NormFz0_new;
%     if NormFz0<=epsilon
%         L1 = 1;
%         NormF = NormFz0; % the final norm of equations
%         break;
%     end
    xik = t*Fz0_newtdk/NormFz0^2;
    z1 = v0-rho*xik*Fz0;
    % compute the next iteration 
    x1 = feval(nprob,z1,2);
    Fx1 = feval(nprob,x1,1);
    NF = NF+1;
    NormFx1 = norm(Fx1);
    if NormFx1<=epsilon
        L1 = 1;
        NormF = NormFx1;
        break;
    end
    
    % update the iteration
   
    x0_old = x0;
    x0 = x1;
    y0_old=v0;
    NormFx0 = NormFx1;
end
if L1==1
    Itr = k;
    Tcpu = toc;
else
    NF = NaN;
    Itr = NaN;
    Tcpu = NaN;
    NormF = NaN;
end
