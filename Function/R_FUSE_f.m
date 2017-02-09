function [X3,A2,time] = R_FUSE_f(H2,M2,E,RE,K,beta,mu,ni,th,...
    ratio,nr,nc,nrr,ncc,Np,Ns,Ne)
%% R-FUSE-TV
tic

dh = zeros(nr,nc); dh(1,1) = 1; dh(1,nc) = -1;
dv = zeros(nr,nc); dv(1,1) = 1; dv(nr,1) = -1;
FDH = fft2(dh); FDHC = conj(FDH);
FDV = fft2(dv); FDVC = conj(FDV);

B=kernel2mat(K,nr,nc);
BF = fft2(B);
BF3 =repmat(BF,[1,1,Ne]);
BF3c=repmat(conj(BF),[1,1,Ne]);
B2Sum=PPlus_s(abs(BF3).^2./(ratio^2),nrr,ncc);

denom=abs(FDH).^2+abs(FDV).^2+1;
IBD_II = 1   ./denom;
IBD_DH = FDHC./denom;
IBD_DV = FDVC./denom;

H3ST=zeros(nr,nc,Ns);
H3ST(1:ratio:end,1:ratio:end,:)=reshape(H2',[nrr,ncc,Ns]);
H2ST=reshape(H3ST,[nr*nc,Ns])';
ETE=E'*E;
[Q,Lambda]=eig(ETE\(RE'*RE+mu*eye(Ne)));
Lambda=reshape(diag(Lambda),[1,1,Ne]);
C4=(ETE*Q)^-1;
Cc=mu*C4;
InvDI=1./(B2Sum+repmat(Lambda,[nrr,ncc,1]));
InvLbd=1./repmat(Lambda,[nr,nc,1]);
C3Cons=fft2(reshape((C4*(E'*H2ST))',[nr,nc,Ne])).*BF3c+...
    fft2(reshape((C4*(RE'*M2))',[nr,nc,Ne]));
C3Cons=C3Cons.*InvLbd;

% initialization
A2=sunsal(E,H2,'lambda',2e-4,'ADDONE','yes','POSITIVITY','yes',...
    'AL_iters',200,'TOL',1e-4,'verbose','no');
[xp,yp,zp]=meshgrid(1:ratio:nc,1:ratio:nr,1:Ne);
[xq,yq,zq]=meshgrid(1:nc,1:nr,1:Ne);
A3i=interp3(xp,yp,zp,reshape(A2',[nrr,ncc,Ne]),xq,yq,zq,'*spline');%'linear'
clear xp yp zp xq yq zq
A2=reshape(A3i,[Np,Ne])';

V2=A2; V3=A2; V4=A2;
G2=zeros(Ne,nr*nc); G3=zeros(Ne,nr*nc); G4=zeros(Ne,nr*nc);
%ci=zeros(ni,1);
%itn=0;

for ii=1:ni
%    itn=itn+1;

    A2=ConvC(V2+G2,IBD_II,nr)+ConvC(V3+G3,IBD_DH,nr)+ConvC(V4+G4,IBD_DV,nr);
    NU2=A2-G2;
    
    VpG=reshape(fft2(reshape((A2-G2)',[nr,nc,Ne])),[nr*nc,Ne])';
    C3bar=C3Cons+reshape((Cc*VpG)',[nr,nc,Ne]).*InvLbd;
    temp=PPlus_s(C3bar/(ratio^2).*BF3,nrr,ncc);
    invQUF=C3bar-repmat(temp.*InvDI,[ratio ratio 1]).*BF3c;
    V2F=Q*reshape(invQUF,[nr*nc,Ne])';
    V2=reshape(real(ifft2(reshape(V2F',[nr,nc,Ne]))),[nr*nc,Ne])';
    
    NU3=ConvC(A2,FDH,nr)-G3;
    NU4=ConvC(A2,FDV,nr)-G4;
    [V3,V4]=vec_soft_col_iso(NU3,NU4,beta/mu);
    G2=-NU2+V2;    % G2 - (X - V2)
    G3=-NU3+V3;    % G3 - (XDh - V3)
    G4=-NU4+V4;    % G4 - (XDv - V4)
    
%    ci(ii)=costc_f(A2,H2,M2,RE,E,K,ratio,nr,nc,nrr,ncc,Ne,alpha,beta);
%    if ii>2 && abs((ci(ii)-ci(ii-1))/ci(ii-1))<th
%        break
%    end
end

% Postprocessing
X3=reshape((E*A2)',[nr,nc,Ns]);
X3(X3>1)=1;
X3(X3<0)=0;

time=toc;