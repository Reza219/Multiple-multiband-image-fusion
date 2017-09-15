function [X3,A2,time] = proposed(H2,M2,p2,E,RE,cE,Kh,Km,Sh,Sm,sp,alpha,...
    mu,ni,ratioh,ratiom,nr,nc,nrh,nch,nrm,ncm,Np,Ns,Ne,Nm)
%% HySure multi lambda (HS+MS+PAN)
tic

proj_simplex_array=@(y)max(bsxfun(@minus,y,max(bsxfun(@rdivide,...
    cumsum(sort(y,1,'descend'),1)-1,(1:size(y,1))'),[],1)),0);

dh=zeros(nr,nc); dh(1,1)=-1; dh(1,nc)=1;
dv=zeros(nr,nc); dv(1,1)=-1; dv(nr,1)=1;
FDh=fft2(dh);
FDv=fft2(dv);

Bh=kernel2mat(Kh,nr,nc);
Bm=kernel2mat(Km,nr,nc);
FBh=fft2(Bh);
FBm=fft2(Bm);

denom=abs(FBh).^2+abs(FBm).^2+abs(FDh).^2+abs(FDv).^2+2;
IBD_I =1   ./denom;

Mh=zeros(nr,nc); %mask, Sh*ShT
Mh(1:ratioh:nr,1:ratioh:nc)=1;
Mh3=repmat(Mh,[1,1,Ne]);
Mh=reshape(Mh3,[nr*nc,Ne])';
H3ST=zeros(nr,nc,Ns);
H3ST(1:ratioh:end,1:ratioh:end,:)=reshape(H2',[nrh,nch,Ns]);
HST=reshape(H3ST,[nr*nc,Ns])';
IEi=(E'*diag(1./Sh)*E+mu*eye(Ne))^-1;
ETHST=E'*diag(1./Sh)*HST;

Mm=zeros(nr,nc); %mask, Sm*SmT
Mm(1:ratiom:nr,1:ratiom:nc)=1;
Mm3=repmat(Mm,[1,1,Ne]);
Mm=reshape(Mm3,[nr*nc,Ne])';
M3ST=zeros(nr,nc,Nm);
M3ST(1:ratiom:end,1:ratiom:end,:)=reshape(M2',[nrm,ncm,Nm]);
MST=reshape(M3ST,[nr*nc,Nm])';
IREi=((RE'*diag(1./Sm)*RE)+mu*eye(Ne))^-1;
ETRTMST=RE'*diag(1./Sm)*MST;

IrTEi=((cE'*cE)/sp+mu*eye(Ne))^-1;
ETrp=(cE'*p2)/sp;

% initialization
U1=sunsal(E,H2,'lambda',2e-4,'ADDONE','yes','POSITIVITY','yes',...
    'AL_iters',200,'TOL',1e-4,'verbose','no');
[xp,yp,zp]=meshgrid(1:ratioh:nc,1:ratioh:nr,1:Ne);
[xq,yq,zq]=meshgrid(1:nc,1:nr,1:Ne);
U1i3=interp3(xp,yp,zp,reshape(U1',[nrh,nch,Ne]),xq,yq,zq,'*spline');%'linear'
clear xp yp zp xq yq zq
U1=reshape(U1i3,[Np,Ne])';

U2=U1; U3=U1; V1=U1; V2=U1; W=U1;
Gu1=zeros(Ne,nr*nc); Gu2=zeros(Ne,nr*nc); Gu3=zeros(Ne,nr*nc);
Gw=zeros(Ne,nr*nc);  Gv1=zeros(Ne,nr*nc); Gv2=zeros(Ne,nr*nc); 

for ii=1:ni
    A2=ConvC(reshape(imfilter(reshape((U1+Gu1)',[nr,nc,Ne]),Kh,'circular'),[Np,Ne])'+...
             reshape(imfilter(reshape((U2+Gu2)',[nr,nc,Ne]),Km,'circular'),[Np,Ne])'+...
             diffchT(V1+Gv1,nr,nc,Np,Ne)+diffcvT(V2+Gv2,nr,nc,Np,Ne)+...
             U3+Gu3+W+Gw,IBD_I,nr);
    Nu1=reshape(imfilter(reshape(A2',[nr,nc,Ne]),Kh,'circular'),[Np,Ne])'-Gu1;
    Nu2=reshape(imfilter(reshape(A2',[nr,nc,Ne]),Km,'circular'),[Np,Ne])'-Gu2;
    Nu3=A2-Gu3;
    Nw =A2-Gw;
    Nv1=diffch(A2,nr,nc,Np,Ne)-Gv1;
    Nv2=diffcv(A2,nr,nc,Np,Ne)-Gv2;
    U1=IEi  *(ETHST  +mu*Nu1.*Mh)+Nu1.*(1-Mh);
    U2=IREi *(ETRTMST+mu*Nu2.*Mm)+Nu2.*(1-Mm);
    U3=IrTEi*(ETrp   +mu*Nu3);
    [V1,V2]=vec_soft_col_iso(Nv1,Nv2,alpha/mu);
    W=proj_simplex_array(Nw);
    Gu1=-Nu1+U1;    % Gu1 - (ABh - U1)
    Gu2=-Nu2+U2;    % Gu2 - (ABm - U2)
    Gu3=-Nu3+U3;    % Gu3 - (A   - U3)
    Gv1=-Nv1+V1;    % Gv1 - (ADh - V1)
    Gv2=-Nv2+V2;    % Gv2 - (ADv - V2)
    Gw =-Nw+W;      % Gu4 - (A   - U4)
end

X3=reshape((E*A2)',[nr,nc,Ns]);
X3(X3>1)=1;
X3(X3<0)=0;

time=toc;