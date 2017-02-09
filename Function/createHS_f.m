function [H3,H2,Kh,nrh,nch,Sh,Noh] = createHS_f(I3,ratioh,sigh,ksh,nr,nc,Ns,Noh,snrh)
%createHS_f creates hyperspectral image

nrh=ceil(nr/ratioh);                      %number of rows
nch=ceil(nc/ratioh);                      %number of columns
Kh=fspecial('gaussian',[ksh,ksh],sigh);   %spatial filter kernel
H3=imfilter(I3,Kh,'circular');
H3=H3(1:ratioh:end,1:ratioh:end,:);
%add noise
%H3=H3+sqrt(mean(H3(:).^2)/10^(snrh/10))*Noh(1:ratioh:end,1:ratioh:end,:);
Sh=zeros(Ns,1);
%Noh=randn(nrh,nch,Ns);
for ii=1:Ns
    Sh(ii,1)=norm(H3(:,:,ii),'fro')^2/(nrh*nch)/10^(snrh/10);
    H3(:,:,ii)=H3(:,:,ii)+sqrt(Sh(ii,1))*Noh(:,:,ii); %1:ratioh:nr,1:ratioh:nc
end
H3(H3>1)=1;   %%%%%%%%%%
H3(H3<0)=0;   %%%%%%%%%%
H2=reshape(H3,[nrh*nch,Ns])';

end

