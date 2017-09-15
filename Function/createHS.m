function [H3,H2,nrh,nch,Sh] = createHS(I3,rh,Kh,nr,nc,Ns,Noh,snrh)
% creates hyperspectral image

nrh=ceil(nr/rh);                      %number of rows
nch=ceil(nc/rh);                      %number of columns
H3=imfilter(I3,Kh,'circular');
H3=H3(1:rh:end,1:rh:end,:);
Sh=zeros(Ns,1);
for ii=1:Ns
    Sh(ii,1)=norm(H3(:,:,ii),'fro')^2/(nrh*nch)/10^(snrh/10);
    H3(:,:,ii)=H3(:,:,ii)+sqrt(Sh(ii,1))*Noh(:,:,ii);
end
H3(H3>1)=1;   %%%%%%%%%%
H3(H3<0)=0;   %%%%%%%%%%
H2=reshape(H3,[nrh*nch,Ns])';

end

