function XD=diffchT(X,nr,nc,N,K)
% circular horizontal gradient

X=reshape(X',[nr,nc,K]);
XD=diff(X,1,2);
XD=cat(2,X(:,1,:)-X(:,end,:),XD);
XD=-reshape(XD,[N,K])';

end