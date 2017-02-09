function [M3,M2,Km,nrm,ncm,Sm,Nom] = createMS_f(I3,ratiom,sigm,ksm,nr,nc,Ns,R,Nm,Nom,snrm)
%createMS_f creates multispectral image

nrm=ceil(nr/ratiom);                      %number of rows
ncm=ceil(nc/ratiom);                      %number of columns
Km=fspecial('gaussian',[ksm,ksm],sigm);   %spatial filter kernel
M3=imfilter(I3,Km,'circular');            % M3*B
M3=M3(1:ratiom:end,1:ratiom:end,:);       % M3*B*D
M2=reshape(M3,[nrm*ncm,Ns])';
M2=R*M2;
% q=max(M2,[],2);
% M2=M2./(q*ones(1,nrm*ncm));
% R =R ./(q*ones(1,Ns));
%add noise
%M2=M2+sqrt(mean(M2(:).^2)/10^(snrm/10))*Nom(1:Nm,1:nrm*ncm);
M3=reshape(M2',[nrm,ncm,Nm]);
Sm=zeros(Nm,1);
%Nom=randn(nrm,ncm,Nm);
for ii=1:Nm
    Sm(ii,1)=norm(M3(:,:,ii),'fro')^2/(nrm*ncm)/10^(snrm/10);
    M3(:,:,ii)=M3(:,:,ii)+sqrt(Sm(ii,1))*Nom(:,:,ii);   % (1:nrm,1:ncm,ii);
end
M3(M3>1)=1;   %%%%%%%%%%
M3(M3<0)=0;   %%%%%%%%%%
M2=reshape(M3,[nrm*ncm,Nm])';

end

