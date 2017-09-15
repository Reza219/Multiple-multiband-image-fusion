function [p3,p2,sp] = createPan(I2,nr,nc,c,nop,snrp)
% creates panchromatic image

p2=c*I2;
sp=mean(p2.^2)/10^(snrp/10);
p2=p2+sqrt(sp)*nop;
p2(p2>1)=1;   %%%%%%%%%%
p2(p2<0)=0;   %%%%%%%%%%
p3=reshape(p2,[nr,nc]);

end

