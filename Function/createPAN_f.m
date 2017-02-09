function [p3,p2,sp,nop] = createPAN_f(I2,nr,nc,c,nop,snrp)
%createMS_f creates panchromatic image

p2=c*I2;
% q=max(p2);
% p2=p2/q;
% c =c /q;
%Nop=randn(1,nr*nc);
%add noise
sp=mean(p2.^2)/10^(snrp/10);
%nop=randn(1,nr*nc);
p2=p2+sqrt(sp)*nop;   % (1:Np)
p2(p2>1)=1;   %%%%%%%%%%
p2(p2<0)=0;   %%%%%%%%%%
p3=reshape(p2,[nr,nc]);

end

