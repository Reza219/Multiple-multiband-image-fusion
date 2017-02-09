function [er,ergas,sam,uiqi,dd] = qual_assess(I3,X3,ratio,Ns,Np)

% relative error
a=squeeze(sum(sum((X3-I3).^2,1),2))/Np;
er=sqrt(sum(a)/Ns)/sqrt(mean(I3(:).^2))*100;

% ERGAS
m=squeeze(sum(sum(I3,1),2))/Np;
ergas=sqrt(sum(a./(m.^2))/Ns)*100/ratio;

% SAM
num=sum(X3.*I3,3);
den=sqrt(sum(X3.^2,3).*sum(I3.^2,3));
sam=sum(acosd(num(:)./den(:)))/Np;

%DD
dd=mean(abs(I3(:)-X3(:)));

%UIQI
q=zeros(Ns,1);
for ii=1:Ns
    q(ii)=img_qi(I3(:,:,ii),X3(:,:,ii),32);
end
uiqi = mean(q);