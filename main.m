clear
addpath('Function')

load('Data\Botswana.mat')
load('Noise\N_Botswana.mat')
alpha=10;                                   %regularization parameter
beta=1.5e-4;                                %regularization parameter
% load('Data\IndianPines.mat')
% load('Noise\N_IndianPines.mat')
% alpha=10;                                   %regularization parameter
% beta=3e-4;                                  %regularization parameter
% load('Data\DCMall.mat')
% load('Noise\N_DCMall.mat')
% alpha=8;                                    %regularization parameter
% beta=1e-3;                                  %regularization parameter
% load('Data\Moffett.mat')
% load('Noise\N_Moffett.mat')
% alpha=15;                                   %regularization parameter
% beta=1e-3;                                  %regularization parameter
% load('Data\KennedySpaceCenter.mat')
% load('Noise\N_KennedySpaceCenter.mat')
% alpha=15;                                   %regularization parameter
% beta=1e-3;                                  %regularization parameter

% paramaters
ni  =200;                                   %number of iterations
snrh=30;                                    %H SNR (dB)
snrm=30;                                    %M SNR (dB)
snrp=40;                                    %p SNR (dB)
mua  =1.5e3;
muo  =.02;                                  %penalty parameter
th  =1e-4;                                  %stopping threshold

%% hyperspectral image
ratioh=4;                           %spatial downsampling factor
sigh  =1.6;
ksh   =7;
[H3,H2,Kh,nrh,nch,Sh]=createHS_f(I3,ratioh,sigh,ksh,nr,nc,Ns,Noh,snrh);

%% spectral responses
R=multibands_Landsat(w);
c=R(6,:);
c=c/sum(c);
R(6,:)=[];   %%%%%%%%%%%
R=R./(sum(R,2)*ones(1,Ns));
RE=R*E;
cE=c*E;

%% multispectral image
ratiom=2;                        %spatial downsampling factor
sigm  =.8;
ksm   =3;
Nm    =7;                        %number of bands
[M3,M2,Km,nrm,ncm,Sm]=createMS_f(I3,ratiom,sigm,ksm,nr,nc,Ns,R,Nm,Nom,snrm);

%% panchromatic
[p3,p2,sp]=createPAN_f(I2,nr,nc,c,nop,snrp);

%% fusion
[X3a,~,tia]=proposed_f(H2,M2,p2,E,RE,cE,Kh,Km,Sh,Sm,sp,alpha,mua,ni,th,...
    ratioh,ratiom,nr,nc,nrh,nch,nrm,ncm,Np,Ns,Ne,Nm);

%% Pan+MS
% BDSD
[M2b,tib]=       BDSD_f(p3,M3,bs,Np,Nm,Km,ratiom);
% MTF-GLP-HPM
[M2m,tim]=MTF_GLP_HPM_f(p3,M3,   Np,Nm,Km,ratiom);

%% MS+HS
[X3hb,~,tihb]=HySure_f(H2,M2b,E,RE,Kh,1,beta,muo,ni,th,ratioh,nr,nc,nrh,nch,Np,Ns,Ne);
[X3fb,~,tifb]=R_FUSE_f(H2,M2b,E,RE,Kh,  beta,muo,ni,th,ratioh,nr,nc,nrh,nch,Np,Ns,Ne);
[X3hm,~,tihm]=HySure_f(H2,M2m,E,RE,Kh,1,beta,muo,ni,th,ratioh,nr,nc,nrh,nch,Np,Ns,Ne);
[X3fm,~,tifm]=R_FUSE_f(H2,M2m,E,RE,Kh,  beta,muo,ni,th,ratioh,nr,nc,nrh,nch,Np,Ns,Ne);

%% results
%Quality indices
[era,ergasa,sama,uiqia,dda]     =qual_assess(I3,X3a ,ratioh,Ns,Np);
[erhb,ergashb,samhb,uiqihb,ddhb]=qual_assess(I3,X3hb,ratioh,Ns,Np);
[erfb,ergasfb,samfb,uiqifb,ddfb]=qual_assess(I3,X3fb,ratioh,Ns,Np);
[erhm,ergashm,samhm,uiqihm,ddhm]=qual_assess(I3,X3hm,ratioh,Ns,Np);
[erfm,ergasfm,samfm,uiqifm,ddfm]=qual_assess(I3,X3fm,ratioh,Ns,Np);

fprintf('proposed\n time = %2.3f\n NRMSE (%%) = %2.3f\n ERGAS = %2.3f\n',...
    tia,era,ergasa)
fprintf(' SAM (degrees) = %2.3f\n UIQI = %2.3f\n DDx10^3 = %2.3f\n \n',...
    sama,uiqia,dda*1e3)

fprintf('BDSD+HySure\n time = %2.3f\n NRMSE (%%) = %2.3f\n ERGAS = %2.3f\n',...
    tihb+tib,erhb,ergashb)
fprintf(' SAM = %2.3f\n UIQI = %2.3f\n DDx10^3 = %2.3f\n \n',...
    samhb,uiqihb,ddhb*1e3)

fprintf('BDSD+RFUSE-TV\n time = %2.3f\n NRMSE (%%) = %2.3f\n ERGAS = %2.3f\n',...
    tifb+tib,erfb,ergasfb)
fprintf(' SAM = %2.3f\n UIQI = %2.3f\n DDx10^3 = %2.3f\n \n',...
    samfb,uiqifb,ddfb*1e3)

fprintf('MTF+HySure\n time = %2.3f\n NRMSE (%%) = %2.3f\n ERGAS = %2.3f\n',...
    tihm+tim,erhm,ergashm)
fprintf(' SAM = %2.3f\n UIQI = %2.3f\n DDx10^3 = %2.3f\n \n',...
    samhm,uiqihm,ddhm*1e3)

fprintf('MTF+RFUSE-TV\n time = %2.3f\n NRMSE (%%) = %2.3f\n ERGAS = %2.3f\n',...
    tifm+tim,erfm,ergasfm)
fprintf(' SAM = %2.3f\n UIQI = %2.3f\n DDx10^3 = %2.3f\n \n',...
    samfm,uiqifm,ddfm*1e3)

%% curves
nrmse_a=zeros(Np,1);
nrmse_hb=zeros(Np,1);
nrmse_fb=zeros(Np,1);
nrmse_hm=zeros(Np,1);
nrmse_fm=zeros(Np,1);
X2a =reshape(X3a, [Np,Ns])';
X2hb=reshape(X3hb,[Np,Ns])';
X2fb=reshape(X3fb,[Np,Ns])';
X2hm=reshape(X3hm,[Np,Ns])';
X2fm=reshape(X3fm,[Np,Ns])';
for ii=1:Np
    nrmse_a(ii) =norm(X2a (:,ii)-I2(:,ii))/norm(I2(:,ii))*100;
    nrmse_hb(ii)=norm(X2hb(:,ii)-I2(:,ii))/norm(I2(:,ii))*100;
    nrmse_fb(ii)=norm(X2fb(:,ii)-I2(:,ii))/norm(I2(:,ii))*100;
    nrmse_hm(ii)=norm(X2hm(:,ii)-I2(:,ii))/norm(I2(:,ii))*100;
    nrmse_fm(ii)=norm(X2fm(:,ii)-I2(:,ii))/norm(I2(:,ii))*100;
end

figure
plot(sort(nrmse_a),'linewidth',3)
hold on
plot(sort(nrmse_hb),'linewidth',3)
plot(sort(nrmse_fb),'linewidth',3)
plot(sort(nrmse_hm),'linewidth',3)
plot(sort(nrmse_fm),'linewidth',3)
legend('proposed','BDSD + HySure','BDSD + RFUSE-TV',...
    'MTF-GLP-HPM + HySure','MTF-GLP-HPM + RFUSE-TV')
xlabel('pixel number')
ylabel('per-pixel NRMSE (%)')

%% images
figure
%subplot(1,5,1)
imshow(p3*gamma)
title('Pan')
axis equal tight
%subplot(1,5,2)
figure
imshow(M3(:,:,4:-1:2)*gamma)
title('MS')
axis equal tight
%subplot(1,5,3)
figure
[~]=imageSPD(H3,w,1.2);
title('HS')
axis equal tight
%subplot(1,5,4)
figure
[~]=imageSPD(I3,w,1.2);
title('original')
axis equal tight
%subplot(1,5,5)
figure
[~]=imageSPD(X3a,w,1.2);
title('fused')
axis equal tight

figure
%subplot(1,2,1)
[~]=imageSPD(X3hb,w,1.2);
title('BDSD+HySure')
axis equal tight
%subplot(1,2,2)
figure
[~]=imageSPD(X3hm,w,1.2);
title('MTF+HySure')
axis equal tight

figure
%subplot(1,2,1)
[~]=imageSPD(X3fb,w,1.2);
title('BDSD+HySureF')
axis equal tight
%subplot(1,2,2)
figure
[~]=imageSPD(X3fm,w,1.2);
title('MTF+HySureF')
axis equal tight