function XD=diffch(X2,nr,nc,Np,Ne)

X3=reshape(X2',[nr,nc,Ne]);
XD=diff(X3,1,2);
XD=cat(2,XD,X3(:,1,:)-X3(:,end,:));
XD=reshape(XD,[Np,Ne])';

end

