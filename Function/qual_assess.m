function [ergas,sam,q2n] = qual_assess(I3,X3,ratio,Ns,Np,Qbs)

% ERGAS
a=squeeze(sum(sum((X3-I3).^2,1),2))/Np;
m=squeeze(sum(sum(I3,1),2))/Np;
ergas=sqrt(sum(a./(m.^2))/Ns)*100/ratio;

% SAM
num=sum(X3.*I3,3);
den=sqrt(sum(X3.^2,3).*sum(I3.^2,3));
sam=nanmean(acosd(num(:)./den(:)));

%Q2n
q2n=q2nf(I3,X3,Qbs);

end