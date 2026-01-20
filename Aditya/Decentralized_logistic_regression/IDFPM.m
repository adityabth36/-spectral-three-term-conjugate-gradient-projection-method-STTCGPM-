% Matlab Model by Jianghua Yin (Jan.,2022, Nanning)
% Copyright (C) 2022 Jian Group
% All Rights Reserved
%
%% the inertial derivative-free projection method (IDFPM) for solving
%% constrained nonlinear pseudo-monotone equations of the form
%   F(x)=0, x\in C,
% where C is a nonempty closed convex set.
%
function [Tcpu,NF,Itr,NormF] = IDFPM(fun,method,Switch,model,x0,x0_old,para)

format long

% start the clock
tic;

%% the number of itrations
% Itr=0;

%% the stopping criterion
epsilon = 1e-6;
epsilon1 = 1e-7;

%% the line search parameters and relaxation factor
k_max = para.Itr_max;   % the maximum number of iterations
gamma = para.gamma;     % the initial guess
sigma = para.sigma;     % the coefficient of line search
tau = para.tau;         % the compression ratio
alpha = para.alpha;     % the coefficient of the inertial step
rho = para.rho;         % the relaxation factor
% alpha_try = 0.1;
mu =10;% para.mu;
lambda_k = 15 ; %para.lambda_k;



fprintf('%s & LSmodel=%d & gamma=%.4f & sigma=%.4f & tau=%.4f & Switch=%.4f & rho=%.4f\n', ...
    method,model,gamma,sigma,tau,Switch,rho);

%% compute the search direction
Fx0 = feval(fun, x0);   % evaluate the function value specified by nprob at x0
NF = 1;
NormFx0 = norm(Fx0);
% x0_old = x0;
L1 = 0;

for k=1:k_max

    if k==1 && NormFx0<=epsilon
        L1 = 1;
        NormF = NormFx0; % the final norm of equations
        break;
    end
    if k==1
        alpha_k = alpha;
    else
        if Switch==1
            alpha_k = min(alpha,1/(norm(x0-x0_old)*k^2));
        elseif Switch==2
            alpha_k = min(alpha,1/(k*norm(x0-x0_old))^2);
        else
            alpha_k = alpha;
        end
    end
    %% compute the inertial step %%
    y0 = x0+alpha_k*(x0-x0_old);
    Fy0 = feval(fun, y0);
    NF = NF+1;
    NormFy0 = norm(Fy0);
    if NormFy0<=epsilon
        L1 = 1;
        NormF = NormFy0;   % the final norm of equations
        break;
    end

    %% compute the initial direction %%
    if k==1
        dk = -Fy0;
    else
        % update the search direction
        switch method
            case 'STTCGPM'
                w0 = Fy0-Fy0_old;
                U_k = mu *(norm(d_k_prev)^2 + NormFy0^2+ abs(d_k_prev'*w0));
                beta_k = NormFy0^2/U_k - ((NormFy0^2 * (Fy0'*d_k_prev))/ U_k^2);
                theta_k= lambda_k+beta_k*((Fy0'*d_k_prev)/(NormFy0^2));
                p_k = Fy0 ;
                dk=-theta_k.*Fy0 + beta_k .*d_k_prev + ((Fy0'*d_k_prev)/(norm(d_k_prev)*norm(p_k)))*p_k;
            case 'FITTCGPM-PRP'
                w0= Fy0-Fy0_old;
                %               betak=(norm(Fk)^2-max(0,norm(Fk)/norm(Fk0)*Fk'*Fk0))/max(norm(Fk0)^2,dk'*yk);
                betak=Fy0'*w0/norm(Fy0_old)^2;
                vk=NormFy0/max(norm(Fy0_old),abs(betak)*norm(dk));
                thetak=-0.055*Fy0'*Fy0_old/norm(Fy0_old)^2;
                dk=-Fy0+0.001*vk*betak*dk+thetak*Fy0_old;
            case'FITTCGPM-DY'%第二篇
                w0 = Fy0-Fy0_old;
                betak=(norm(Fy0))^2/(dk'*w0);
                nuk=norm(Fy0)/max(norm(Fy0_old),(abs(betak))*norm(dk));
                thetak=-0.055*Fy0'*Fy0_old/norm(Fy0_old)^2;  % 原始的分母是 norm(Fk0)^2
                dk=-Fy0+0.001*nuk*betak*dk+thetak*Fy0_old;
            case'FITTCGPM-OUR'
                w0 = Fy0-Fy0_old;
                U_k = mu *(norm(d_k_prev)^2 + NormFy0^2+ abs(d_k_prev'*w0));
                betak = NormFy0^2/U_k - ((NormFy0^2 * (Fy0'*d_k_prev))/ U_k^2);
                 nuk=norm(Fy0)/max(norm(Fy0_old),(abs(betak))*norm(dk));
                thetak=-0.055*Fy0'*Fy0_old/norm(Fy0_old)^2;  % 原始的分母是 norm(Fk0)^2
                dk=-Fy0+0.001*nuk*betak*dk+thetak*Fy0_old;
             case 'IMSMNE'
                t1=0.3; t2=10^-3; D=10; r0=3; tauk=6.2 ; 
                w0 = Fy0-Fy0_old;
                sk=y0-y0_old;
                h_k = w0 + D*(norm(Fy0_old)^r0) *sk+ max(0,((-sk'*w0)/(sk'*sk)))*sk;
                k_a = max(tauk*norm(d_k_prev)*norm(h_k), d_k_prev'*h_k);
                betak= (Fy0'*h_k)/(k_a) - ((norm(h_k)^2 * (Fy0'*d_k_prev))/(k_a^2))  ;
                theta = ((sk'*h_k)/(sk'*sk));
                theta_k =  min(max(theta, t1), t2);
                dk = - theta_k * Fy0 + betak * d_k_prev ;     
            otherwise
                disp('Input error! Please check the input method');
        end
    end
    Normdk = norm(dk);
    if Normdk<epsilon1
        L1 = 1;
        NormF = NormFy0;
        break;
    end
    Normdk2 = Normdk^2;
    Fy0_old = Fy0;

    %%% Start Armijo-type line search  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % model=1 means -F(zk)'*dk ?? sigma*tk*norm(dk)^2
    % model=2 means -F(zk)'*dk ?? sigma*tk*norm(F(zk))*norm(dk)^2
    % model=3 means -F(zk)'*dk ?? sigma*tk*norm(F(zk))/(1+norm(F(zk)))*norm(dk)^2
    % model=4 means -F(zk)'*dk ?? sigma*tk*max(lambda,min(nu,norm(Fz_new,2)))*norm(dk)^2
    if model==1
        t = gamma;
        z0_new = y0+t*dk;
        Fz0_new = feval(fun, z0_new);
        NF = NF+1;
        Fz0_newtdk = -Fz0_new'*dk;
        % check the Armijo-type line search condition
        while Fz0_newtdk < sigma*t*Normdk2 %&& t>10^-6
            % the Armijo-type line search condition violated
            t = t*tau;
            z0_new = y0+t*dk;
            Fz0_new = feval(fun, z0_new);
            NF = NF+1;
            Fz0_newtdk = -Fz0_new'*dk;
        end %%% End Armijo-type line search %%%
        NormFz0_new = norm(Fz0_new);
    elseif model==2
        t = gamma;
        z0_new = y0+t*dk;
        Fz0_new = feval(fun, z0_new);
        NF = NF+1;
        NormFz0_new = norm(Fz0_new);
        Fz0_newtdk = -Fz0_new'*dk;
        % check the Armijo-type line search condition
        while Fz0_newtdk < sigma*t*NormFz0_new*Normdk2 %&& t>10^-6
            % the Armijo-type line search condition violated
            t = t*tau;
            z0_new = y0+t*dk;
            Fz0_new = feval(fun, z0_new);
            NF = NF+1;
            NormFz0_new = norm(Fz0_new);
            Fz0_newtdk = -Fz0_new'*dk;
        end %%% End Armijo-type line search %%%
    elseif model==3
        t = gamma;
        z0_new = y0+t*dk;
        Fz0_new = feval(fun, z0_new);
        NF = NF+1;
        NormFz0_new = norm(Fz0_new);
        Fz0_newtdk = -Fz0_new'*dk;
        % check the Armijo-type line search condition
        while Fz0_newtdk < sigma*t*NormFz0_new/(1+NormFz0_new)*Normdk2 && t>10^-6
            % the Armijo-type line search condition violated
            t = t*tau;
            z0_new = y0+t*dk;
            Fz0_new = feval(fun, z0_new);
            NF = NF+1;
            NormFz0_new = norm(Fz0_new);
            Fz0_newtdk = -Fz0_new'*dk;
        end %%% End Armijo-type line search %%%
    else
        t = gamma;
        z0_new = y0+t*dk;
        Fz0_new = feval(fun, z0_new);
        NF = NF+1;
        NormFz0_new = norm(Fz0_new);
        Fz0_newtdk = -Fz0_new'*dk;
        % check the Armijo-type line search condition
        while Fz0_newtdk < sigma*t*max(0.001,min(0.8,NormFz0_new))*Normdk2 && t>10^-6
            % the Armijo-type line search condition violated
            t = t*tau;
            z0_new = y0+t*dk;
            Fz0_new = feval(fun, z0_new);
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
    % compute the next iteration
    x1 = y0-rho*xik*Fz0;
    Fx1 = feval(fun, x1);
    NF = NF+1;
    NormFx1 = norm(Fx1);
    if NormFx1<=epsilon
        L1 = 1;
        NormF = NormFx1;
        break;
    end

    % update the iteration
    d_k_prev =dk;
    x0_old = x0;
    y0_old = y0;
    x0 = x1;
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