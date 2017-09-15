function [M2f,time] = MTF_GLP_HPM(p3,M3,Np,Nm,Km,ratio)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description:
%           MTF_GLP_HPM fuses the upsampled MultiSpectral (MS) and PANchromatic (PAN) images by
%           exploiting the Modulation Transfer Function - Generalized Laplacian Pyramid (MTF-GLP)
%           with High Pass Modulation (HPM) injection model algorithm.
%
% Inputs:
%           p3:              PAN image;
%           M3:              MS image;
%           ratio:           Scale ratio between MS and PAN.
%
% Outputs:
%           M2f:             pansharpened MS image.
%
% References:
%           [Aiazzi03]       B. Aiazzi, L. Alparone, S. Baronti, A. Garzelli, and M. Selva, “An MTF-based spectral distortion minimizing model for Pan-sharpening
%                            of very high resolution multispectral images of urban areas,” in Proceedings of URBAN 2003: 2nd GRSS/ISPRS Joint Workshop on
%                            Remote Sensing and Data Fusion over Urban Areas, 2003, pp. 90–94.
%           [Aiazzi06]       B. Aiazzi, L. Alparone, S. Baronti, A. Garzelli, and M. Selva, “MTF-tailored multiscale fusion of high-resolution MS and Pan imagery,”
%                            Photogrammetric Engineering and Remote Sensing, vol. 72, no. 5, pp. 591–596, May 2006.
%           [Vivone14a]      G. Vivone, R. Restaino, M. Dalla Mura, G. Licciardi, and J. Chanussot, “Contrast and error-based fusion schemes for multispectral
%                            image pansharpening,” IEEE Geoscience and Remote Sensing Letters, vol. 11, no. 5, pp. 930–934, May 2014.
%           [Vivone14b]      G. Vivone, L. Alparone, J. Chanussot, M. Dalla Mura, A. Garzelli, G. Licciardi, R. Restaino, and L. Wald, “A Critical Comparison Among Pansharpening Algorithms”,
%                            IEEE Transaction on Geoscience and Remote Sensing, 2014. (Accepted)
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
p3l=imfilter(p3h,Km,'circular');
t  =p3l(1:ratio:end,1:ratio:end,:);
p3l=interp23tap(t,ratio);

M3f=M3u.*(p3h./(p3l+eps));
M2f=reshape(M3f,[Np,Nm])';

time=toc;
end