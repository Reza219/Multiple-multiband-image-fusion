function [M2f,time] = BDSD_f(p3,M3,S,Np,Nm,Km,ratio)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description: 
%           BDSD fuses the upsampled MultiSpectral (MS) and PANchromatic (PAN) images by 
%           exploiting the Band-Dependent Spatial-Detail (BDSD) algorithm. 
% 
% Interface:
%           I_Fus_BDSD = BDSD(I_MS,I_PAN,ratio,S,sensor)
%
% Inputs:
%           I_MS:           MS image upsampled at PAN scale;
%           I_PAN:          PAN image;
%           ratio:          Scale ratio between MS and PAN. Pre-condition: Integer value;
%           S:              Local estimation on SxS distinct blocks (typically 128x128); 
%           sensor:         String for type of sensor (e.g. 'WV2', 'IKONOS').
%
% Outputs:
%           I_Fus_BDSD:     BDSD pansharpened image.
% 
% References:
%           [Garzelli08]    A. Garzelli, F. Nencini, and L. Capobianco,
%                           “Optimal MMSE pan sharpening of very high
%                           resolution multispectral images,” IEEE
%                           Transactions on Geoscience and Remote Sensing,
%                           vol. 46, no. 1, pp. 228–236, January 2008.
%           [Vivone14]      G. Vivone, L. Alparone, J. Chanussot,
%                           M. Dalla Mura, A. Garzelli, G. Licciardi,
%                           R. Restaino, and L. Wald, “A Critical
%                           Comparison Among Pansharpening Algorithms”, 
%                           IEEE Transaction on Geoscience and Remote
%                           Sensing, 2014. (Accepted)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic
M3u=interp23tap(M3,ratio);

%%%
% Control of input parameters and initialization
%%%
if S>1
    if(rem(S,2) && S >1)
        fprintf(1,'\n\n ');
        error('block size for local estimation must be even')
    end

    if(rem(S,ratio))
        fprintf(1,'\n\n ');
        error('block size must be multiple of ratio')
    end

    [nr,nc] = size(p3);

    if(rem(nr,S)||rem(nc,S))
        fprintf(1,'\n\n ');
        error('x and y dims of pan must be multiple of the block size')
    end
end

%%%
% Reduced resolution
%%%

% p3l=MTF_PAN(p3,sensor,ratio);
p3l=imfilter(p3,Km,'circular');
p3l=p3l(1:ratio:end,1:ratio:end);
% M3l=MTF(M3,sensor,[],ratio);
M3l=zeros(size(M3));
for ii=1:Nm
    M3l(:,:,ii)=imfilter(M3(:,:,ii),Km,'circular');
end

%%%
% Parameter estimation at reduced resolution
%%%
in3 = cat(3,M3l,M3,p3l);
fun_eg = @(bs) estimate_gamma_cube(bs.data,S,ratio);
gamma = blockproc(in3,[S/ratio S/ratio],fun_eg);

%%%
% Fusion
%%%
in3 = cat(3,M3u,p3,gamma);
fun_Hi = @(bs) compH_inject(bs.data,S);
M3f = blockproc(in3,[S,S],fun_Hi);

M2f=reshape(M3f,[Np,Nm])';
time=toc;
end

%%%_______________________________________________________________
%%%
function gamma = estimate_gamma_cube(in3,S,ratio)
Nb = (size(in3,3)-1)/2;
hs_LP_d = in3(:,:,1:Nb);
hs_orig = in3(:,:,Nb+1:2*Nb);
pan_LP_d = in3(:,:,2*Nb+1);
% Compute Hd
Hd = zeros(S*S/ratio/ratio,Nb+1);
for k=1:Nb
    b = hs_LP_d(:,:,k);
    Hd(:,k) = b(:);
end
Hd(:,Nb+1) = pan_LP_d(:);
% Estimate gamma
B = (Hd'*Hd)\Hd';
gamma = zeros(Nb+1,Nb);
for k=1:Nb
    b = hs_orig(:,:,k);
    bd = hs_LP_d(:,:,k);
    gamma(:,k) = B *(b(:)-bd(:));
end
gamma = padarray(gamma,[S-Nb-1 S-Nb],0,'post');
end

%%%_______________________________________________________________
%%%
function ms_en = compH_inject(in3,S)
Nb = size(in3,3)-2;
hs = in3(:,:,1:Nb);
pan = in3(:,:,Nb+1);
gamma = in3(:,:,Nb+2); 
% Compute H
[N,M,Nb] = size(hs);
H = zeros(S*S,Nb+1);
for k=1:Nb
    b = hs(:,:,k);
    H(:,k) = b(:);
end
H(:,Nb+1) = pan(:);
% Inject
g = gamma(1:Nb+1,1:Nb);
ms_en = zeros(N,M,Nb);
for k=1:Nb
    b = hs(:,:,k);
    b_en = b(:) + H * g(:,k);
    ms_en(:,:,k) = reshape(b_en,N,M);
end
end