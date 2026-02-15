% This is a demo of DFPM for solving \ell_1-norm minimization
% problem
% arg min_x = 0.5*|| A x - b ||_2^2 + lambda || x ||_1
% using the algorithm modified three-term conjugate gradient projection method, described in the following paper



clc
clear all
close all
clf
randn('seed',1); % 1 for the experiments in the paper
rand('seed',1); % 1 for the experiments in the paper
addpath Images   %It adds the folder named Images to the MATLAB search path
% the test images set
A={'fruits.bmp' 'baboon.bmp' 'brain.bmp' 'cameraman.png' 'lena.png' 'peppers.png' 'circles.png'...
    'chart.png' 'barbara.png' 'man.bmp' 'shape.jpg' 'housergb.png' 'house.png'};
% f = double(imread('cameraman.png'));
% f = double(imread('lena.png'));
% f = double(imread('barbara.png'));
% f = double(imread('man.bmp'));
f=double(imread(A{8})); 
[m,n] = size(f);
scrsz = get(0,'ScreenSize');
figure(1)
%set(1,'Position',[10 scrsz(4)*0.05 scrsz(3)/4 0.85*scrsz(4)])
% subplot(1,2,1)
imagesc(f)
colormap(gray(255))
axis off
axis equal
% title('Original image','FontName','Times','FontSize',22)
% title(sprintf('Original: %4d  %4d',m,n),'FontSize',22)
% create observation operator; in this case 
% it will be a blur function composed with an
% inverse weavelet transform
disp('Creating observation operator...');

middle = n/2 + 1;

% uncomment the following lines for Experiment 1 (see paper)
% sigma = 0.56;
% h = zeros(size(f));
% for i=-4:4
%    for j=-4:4
%       h(i+middle,j+middle)= 1; 
%    end
% end


% uncomment the following lines for Experiment 2 (see paper)
% sigma = sqrt(2);
sigma = 0.001;
h = zeros(size(f)); % h is the same size as f and all zeros.
for i=-4:4
   for j=-4:4
      h(i+middle,j+middle)= (1/(1+i*i+j*j));
   end
end

% uncomment the following lines for Experiment 3 (see paper)
% sigma = sqrt(8);
% h = zeros(size(f));
% for i=-4:4
%    for j=-4:4
%       h(i+middle,j+middle)= (1/(1+i*i+j*j));
%    end
% end

% % center and normalize the blur
h = fftshift(h);  % ½«ÆµÆ×Í¼ÒÆÎ»£¬µÍÆµÒÆÖÁÆµÆ×Í¼ÖÐÐÄ
% fftshift is useful for visualizing the Fourier transform with the zero-frequency component in the middle of the spectrum.
h = h/sum(h(:));  % sum(h(:)): sum all elements of h.

% definde the function handles that compute 
% the blur and the conjugate blur.
R = @(x) real(ifft2(fft2(h).*fft2(x))); % fft2(X) returns the two-dimensional Fourier transform of matrix X.
% ifft2(F) returns the two-dimensional inverse Fourier transform of matrix F
RT = @(x) real(ifft2(conj(fft2(h)).*fft2(x)));

% define the function handles that compute 
% the products by W (inverse DWT) and W' (DWT)
wav = daubcqf(2);
W = @(x) midwt(x,wav,3);
WT = @(x) mdwt(x,wav,3);

%Finally define the function handles that compute 
% the products by A = RW  and A' =W'*R' 
A = @(x) R(W(x));
AT = @(x) WT(RT(x));

fid=fopen('mytext.txt','w');
% generate noisy blurred observations
b = R(f) + sigma*randn(size(f));
SNR_original = 20*log10(norm(f,'fro')/norm(f-b,'fro'));
figure(2)
% subplot(1,2,2)
imagesc(b)
colormap(gray(255))
axis off
axis equal

ITR_max = 5000;
% set parameters for IMSMNE
para4.Itr_max = ITR_max;
para4.gamma = 1;         % the initial guess
para4.sigma = 0.001;         % the coefficient of line search 
para4.tau = 0.36;         % the compression ratio
para4.alpha = 0.35;       % the coefficient of inertial step
para4.rho = 1.87;         % the relaxation factor 



% parameters for ISTCP
para5.Itr_max = ITR_max;
para5.gamma = 1;         % the initial guess
para5.sigma = 0.0001;      % the coefficient of line search 
para5.tau = 0.7;         % the compression ratio
para5.alpha = 0.8;       % the coefficient of inertial step
para5.rho = 1;      % the relaxation factor 


% parameters for STTCGPM
% para6.Itr_max = ITR_max;
% para6.gamma = 0.4;         % the initial guess
% para6.sigma = 0.0001;      % the coefficient of line search 
% para6.tau = 0.3;         % the compression ratio
% para6.alpha = 0.5;      % the coefficient of inertial step
% para6.rho = 1.9;      % the relaxation factor 

% para6.Itr_max = ITR_max;
% para6.gamma = 0.72;         % the initial guess
% para6.sigma = 0.01;      % the coefficient of line search 
% para6.tau = 0.6;         % the compression ratio
% para6.alpha = 0.29;      % the coefficient of inertial step
% para6.rho = 1.8882;      % the relaxation factor 
para6.Itr_max = ITR_max;
para6.gamma = 0.4;         % the initial guess
para6.sigma = 0.01;      % the coefficient of line search 
para6.tau = 0.3;         % the compression ratio
para6.alpha = 0.3;      % the coefficient of inertial step
para6.rho = 1.9;      % the relaxation factor 


%FITTCGPM-PRP
para7.Itr_max = ITR_max;
para7.gamma = 0.72;         % the initial guess0.6
para7.sigma = 0.01;      % the coefficient of line search 
para7.tau = 0.6;         % the compression ratio
para7.alpha = 0.29;      % the coefficient of inertial step
para7.rho = 1.8882;      % the relaxation factor 


% regularization parameter
varrho = 0.035; %0.35 .0525 .0635 0.035

disp('Starting STTCGPM')
[x1,SNR1,Tcpu1,NF1,NormF1] = IRDFPM(A,AT,W,b,varrho,f,'STTCGPM',2,para6);
T1 = Tcpu1(end);

disp('Starting FITTCGPM-PRP')
[x2,SNR2,Tcpu2,NF2,NormF2] = IRDFPM(A,AT,W,b,varrho,f,'FITTCGPM-PRP',1,para7);
T2 = Tcpu2(end);
% USDGPM_mses = mses2(end)

disp('Starting FITTCGPM-DY')
[x3,SNR3,Tcpu3,NF3,NormF3] = IRDFPM(A,AT,W,b,varrho,f,'FITTCGPM-DY',1,para7);
T3 = Tcpu3(end);
% USDGPM_mses = mses2(end)

disp('Starting IMSMNE')
[x4,SNR4,Tcpu4,NF4,NormF4] = IRDFPM(A,AT,W,b,varrho,f,'IMSMNE',4,para4);
T4 = Tcpu4(end);
% IRSDGPM_mses = mses1(end)




% ================= Plotting results =================
%%  semilogy
fprintf(fid,'%.2f & %.2f/%.2f & %.2f/%.2f & %.2f/%.2f & %.2f/%.2f\n', ... 
        SNR_original,T1(end),SNR1(end),T2(end),SNR2(end),T3(end),SNR3(end),T4(end),SNR4(end));
fclose(fid);

figure(3)
% subplot(1,5,4)
imagesc(W(x1))
colormap(gray)
axis off
axis equal
% title(['YJJ HCGP: ', num2str(time_YJJ_HCGP),'s',num2str(SNR_YJJ_HCGP),'dB'], 'FontName','Times','FontSize',22)

figure(4)
% subplot(1,5,4)
imagesc(W(x2))
colormap(gray)
axis off
axis equal
% title(['FITTCGPM_DY: ', num2str(time_SGP),'s',num2str(SNR_SGP),'dB'], 'FontName','Times','FontSize',22)

figure(5)
% subplot(1,5,4)
imagesc(W(x3))
colormap(gray)
axis off
axis equal
% title(['HCGP: ', num2str(time_HCGP),'s',num2str(SNR_HCGP),'dB'], 'FontName','Times','FontSize',22)

figure(6)
% subplot(1,5,4)
imagesc(W(x4))
colormap(gray)
axis off
axis equal

% figure(7)
% % subplot(1,5,4)
% imagesc(W(x5))
% colormap(gray)
% axis off
% axis equal

% % scrsz = get(0,'ScreenSize');
% % set(7,'Position',[10 scrsz(4)*0.1 0.9*scrsz(3)/2 3*scrsz(4)/4])
% subplot(5,1,5)
% plot(x4(:),'LineWidth',1.1)
% top = max(x4(:));
% bottom = min(x4(:));
% v = [0 n+1 bottom-0.05*(top-bottom)  top+0.05*((top-bottom))];
% set(gca,'FontName','Times')
% set(gca,'FontSize',14)
% title(sprintf('TTCGPM (MSE = %5.2e, Itr=%g, Tcpu=%4.2fs)',mses4(end),length(mses4),Tcpu4(end)))
% axis(v)
% % 
% % 
