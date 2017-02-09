function [X3,A2,time] = HySure_f(H2,M2,E,RE,K,alpha,beta,mu,ni,th,ratio,...
    nr,nc,nrr,ncc,Np,Ns,Ne)
%% HySure
tic

dh = zeros(nr,nc); dh(1,1) = 1; dh(1,nc) = -1;
dv = zeros(nr,nc); dv(1,1) = 1; dv(nr,1) = -1;
FDH = fft2(dh); FDHC = conj(FDH);
FDV = fft2(dv); FDVC = conj(FDV);

B=kernel2mat(K,nr,nc);
FB = fft2(B); FBC = conj(FB);

denom=abs(FB.^2)+abs(FDH).^2+abs(FDV).^2+1;
IBD_B  = FBC ./denom;
IBD_II = 1   ./denom;
IBD_DH = FDHC./denom;
IBD_DV = FDVC./denom;

mask = zeros(nr, nc);
mask(1:ratio:nr, 1:ratio:nc) = 1;
maskim = repmat(mask, [1, 1, Ne]);
mask = reshape(maskim,[nr*nc,Ne])';

H3ST=zeros(nr,nc,Ns);
H3ST(1:ratio:end,1:ratio:end,:)=reshape(H2',[nrr,ncc,Ns]);
H2ST = reshape(H3ST,[nr*nc,Ns])';
IE = E'*E+mu*eye(Ne);
yyh = E'*H2ST;
IRE = alpha*(RE'*RE)+mu*eye(Ne);
yym = RE'*M2;

% initialization
V1=sunsal(E,H2,'lambda',2e-4,'ADDONE','yes','POSITIVITY','yes',...
    'AL_iters',200,'TOL',1e-4,'verbose','no');
[xp,yp,zp]=meshgrid(1:ratio:nc,1:ratio:nr,1:Ne);
[xq,yq,zq]=meshgrid(1:nc,1:nr,1:Ne);
V1i3=interp3(xp,yp,zp,reshape(V1',[nrr,ncc,Ne]),xq,yq,zq,'*spline');%'linear'
clear xp yp zp xq yq zq
V1=reshape(V1i3,[Np,Ne])';

V2=V1; V3=V1; V4=V1;
G1=zeros(Ne,nr*nc); G2=zeros(Ne,nr*nc); G3=zeros(Ne,nr*nc); G4=zeros(Ne,nr*nc);
%ci=zeros(ni,1);
%itn=0;

for ii=1:ni
%    itn=itn+1;
    
    A2 = ConvC(V1+G1, IBD_B, nr) + ConvC(V2+G2, IBD_II, nr) + ...
        ConvC(V3+G3, IBD_DH, nr) +  ConvC(V4+G4, IBD_DV, nr);
    NU1 =  ConvC(A2, FB, nr) - G1;
    V1 = IE\(yyh + mu*NU1).*mask + NU1.*(1-mask);
    NU2 =  A2 - G2;
    V2 = IRE\(alpha*yym + mu*NU2);
    NU3 =  ConvC(A2, FDH, nr) - G3;
    NU4 =  ConvC(A2, FDV, nr) - G4;
    [V3,V4] = vec_soft_col_iso(NU3,NU4,beta/mu);
    G1 = -NU1 + V1;    % D1 - (XB - V1)
    G2 = -NU2 + V2;    % D2 - (XDH - V2)
    G3 = -NU3 + V3;    % D3 - (XDV - V3)
    G4 = -NU4 + V4;    % D3 - (XDV - V3)
    
%    ci(ii)=costc_f(A2,H2,M2,RE,E,K,ratio,nr,nc,nrr,ncc,Ne,alpha,beta);
%    if ii>2 && abs((ci(ii)-ci(ii-1))/ci(ii-1))<th
%        break
%    end
end

X3=reshape((E*A2)',[nr,nc,Ns]);
X3(X3>1)=1;
X3(X3<0)=0;

time=toc;