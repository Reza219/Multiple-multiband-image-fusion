function XD=diffcv_f(X2,nr,nc,Np,Ne)

X3=reshape(X2',[nr,nc,Ne]);
XD=diff(X3);
XD=cat(1,XD,X3(1,:,:)-X3(end,:,:));
XD=reshape(XD,[Np,Ne])';

end

