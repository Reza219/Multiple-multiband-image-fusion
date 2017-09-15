clear
addpath('function')

load('Data\Botswana.mat')
load('Noise\N_Botswana.mat')
beta=6e-5;                                  % regularization parameter
ps=14:36;
% load('Data\IndianPines.mat')
% load('Noise\N_IndianPines.mat')
% beta=7e-5;                                  % regularization parameter
% ps=10:30;
% load('Data\DCMall.mat')
% load('Noise\N_DCMall.mat')
% beta=2.5e-4;                                % regularization parameter
% ps=24:53;
% load('Data\Moffett.mat')
% load('Noise\N_Moffett.mat')
% beta=8e-4;                                  % regularization parameter
% ps=14:36;
% load('Data\KennedySpaceCenter.mat')
% load('Noise\N_KennedySpaceCenter.mat')
% beta=8e-4;                                  % regularization parameter
% ps=14:36;

% paramaters
ni  =200;                                   % number of iterations
snrh=30;                                    % HS SNR (dB)
snrm=30;                                    % MS SNR (dB)
snrp=40;                                    % pan SNR (dB)
mu  =.02;                                   % penalty parameter

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% spectral responses
R=multibands_Landsat(w);
c=R(6,:);
c=c/sum(c);
cE=c*E;

%% MTFs (blur kernels)
Km =fspecial('gaussian',[3,3],.8);       % MS spatial filter kernel
rh =4;                                   % HS spatial downsampling factor
Kh =conv2(Km,Km);                        % HS spatial filter kernel

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% synthesize images
[~,p2,~]         =createPan(I2,nr,nc,c,nop,snrp);               % Pan
[~,H2,nrh,nch,Sh]=createHS (I3,rh,Kh,nr,nc,Ns,Noh,snrh);        % HS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Pan+HS fusion
[X3,~,time]=HySure(H2,p2,E,cE,Kh,1,beta,mu,ni,rh,nr,nc,nrh,nch,Np,Ns,Ne);
%[X3,~,time]=RFUSE_TV(H2,p2,E,cE,Kh,beta,mu,ni,rh,nr,nc,nrh,nch,Np,Ns,Ne);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% results
[ergas,sam,q2n]=qual_assess(I3(:,:,ps),X3(:,:,ps),rh,length(ps),Np,32);
fprintf('\n spectrum of pan\n ERGAS = %2.3f\n SAM   = %2.3f\n Q2n   = %2.3f',...
    ergas,sam,q2n)

[ergas,sam,q2n]=qual_assess(I3,X3,rh,Ns,Np,32);
fprintf('\n entire spectrum\n ERGAS = %2.3f\n SAM   = %2.3f\n Q2n   = %2.3f\n time  = %2.3f\n'...
    ,ergas,sam,q2n,time)

nrmse=zeros(Np,1);
X2p=reshape(X3(:,:,ps),[Np,length(ps)])';
I2p=reshape(I3(:,:,ps),[Np,length(ps)])';
for ii=1:Np
    nrmse(ii)=norm(X2p(:,ii)-I2p(:,ii))/norm(I2p(:,ii))*100;
end
figure, plot(sort(nrmse)), title('spectrum of pan')
xlabel('pixel number'), ylabel('per-pixel NRMSE (%)')

nrmse=zeros(Np,1);
X2=reshape(X3,[Np,Ns])';
for ii=1:Np
    nrmse(ii)=norm(X2(:,ii)-I2(:,ii))/norm(I2(:,ii))*100;
end
figure, plot(sort(nrmse)), title('entire spectrum')
xlabel('pixel number'), ylabel('per-pixel NRMSE (%)')