function [M2f,time] = MTF_GLP_HPM_f(p3,M3,Np,Nm,Km,ratio)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description:
%           MTF_GLP_HPM fuses the upsampled MultiSpectral (MS) and PANchromatic (PAN) images by
%           exploiting the Modulation Transfer Function - Generalized Laplacian Pyramid (MTF-GLP)
%           with High Pass Modulation (HPM) injection model algorithm.
%
% Interface:
%           I_Fus_MTF_GLP_HPM = MTF_GLP_HPM(I_PAN,I_MS,sensor,tag,ratio)
%
% Inputs:
%           I_PAN:              PAN image;
%           I_MS:               MS image upsampled at PAN scale;
%           sensor:             String for type of sensor (e.g. 'WV2','IKONOS');
%           tag:                Image tag. Often equal to the field sensor.
%                               It makes sense when sensor is 'none'. It indicates the band number
%                               in the latter case;
%           ratio:              Scale ratio between MS and PAN. Pre-condition: Integer value.
%
% Outputs:
%           I_Fus_MTF_GLP_HPM:  MTF_GLP_HPM pansharpened image.
%
% References:
%           [Aiazzi03]          B. Aiazzi, L. Alparone, S. Baronti, A. Garzelli, and M. Selva, “An MTF-based spectral distortion minimizing model for Pan-sharpening
%                               of very high resolution multispectral images of urban areas,” in Proceedings of URBAN 2003: 2nd GRSS/ISPRS Joint Workshop on
%                               Remote Sensing and Data Fusion over Urban Areas, 2003, pp. 90–94.
%           [Aiazzi06]          B. Aiazzi, L. Alparone, S. Baronti, A. Garzelli, and M. Selva, “MTF-tailored multiscale fusion of high-resolution MS and Pan imagery,”
%                               Photogrammetric Engineering and Remote Sensing, vol. 72, no. 5, pp. 591–596, May 2006.
%           [Vivone14a]         G. Vivone, R. Restaino, M. Dalla Mura, G. Licciardi, and J. Chanussot, “Contrast and error-based fusion schemes for multispectral
%                               image pansharpening,” IEEE Geoscience and Remote Sensing Letters, vol. 11, no. 5, pp. 930–934, May 2014.
%           [Vivone14b]         G. Vivone, L. Alparone, J. Chanussot, M. Dalla Mura, A. Garzelli, G. Licciardi, R. Restaino, and L. Wald, “A Critical Comparison Among Pansharpening Algorithms”,
%                               IEEE Transaction on Geoscience and Remote Sensing, 2014. (Accepted)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic
M3u=interp23tap(M3,ratio);

%%% Equalization
p3h=repmat(p3,[1,1,Nm]);
for ii=1:Nm
    p3h(:,:,ii)=(p3h(:,:,ii)-mean2(p3h(:,:,ii))).*...
        (std2(M3u(:,:,ii))./std2(p3h(:,:,ii)))+mean2(M3u(:,:,ii));
end

%%% MTF
%N  =41;
p3l=zeros(size(M3u));

for ii=1:Nm
    %alpha      =sqrt((N/ratio/2)^2/(-2*log(GNyq(ii))));
    %H          =fspecial('gaussian',N,alpha);
    %H          =H./max(H(:));
    %h          =fwind1(H,kaiser(N));
    %p3l(:,:,ii)=imfilter(M3f(:,:,ii),real(h),'replicate');%'symmetric''circular'
    p3l(:,:,ii)=imfilter(p3h(:,:,ii),Km,'circular');%'symmetric''circular'
    %t          =imresize(p3l(:,:,ii),1/ratio);%,'nearest');
    t          =p3l(1:ratio:end,1:ratio:end,ii);
    %t          =p3l(2:ratio:end,2:ratio:end,ii);
    p3l(:,:,ii)=interp23tap(t,ratio);
end

M3f=M3u.*(p3h./(p3l+eps));
M2f=reshape(M3f,[Np,Nm])';

time=toc;
end