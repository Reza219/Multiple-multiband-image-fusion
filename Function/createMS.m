function [M3,M2,nrm,ncm,Sm] = createMS(I3,rm,Km,nr,nc,Ns,R,Nm,Nom,snrm)
% creates multispectral image

nrm=ceil(nr/rm);                      %number of rows
ncm=ceil(nc/rm);                      %number of columns
M3=imfilter(I3,Km,'circular');
M3=M3(1:rm:end,1:rm:end,:);
M2=reshape(M3,[nrm*ncm,Ns])';
M2=R*M2;
M3=reshape(M2',[nrm,ncm,Nm]);
Sm=zeros(Nm,1);
for ii=1:Nm
    Sm(ii,1)=norm(M3(:,:,ii),'fro')^2/(nrm*ncm)/10^(snrm/10);
    M3(:,:,ii)=M3(:,:,ii)+sqrt(Sm(ii,1))*Nom(:,:,ii);
end
M3(M3>1)=1;   %%%%%%%%%%
M3(M3<0)=0;   %%%%%%%%%%
M2=reshape(M3,[nrm*ncm,Nm])';

end

