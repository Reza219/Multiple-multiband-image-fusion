%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           [Q2n_index,Q2n_index_map] = q2n(I3,X3,Q_block_size)
%
% Inputs:
%           I3:              Ground-Truth image;
%           X3:              Fused Image;
%           bs:              Block size of the Q-index locally applied;
%
% Outputs:
%           Q2n_index:       Q2n index;
%           Q2n_map:         Map of Q2n values.
%
% References:
%           [Garzelli09]     A. Garzelli and F. Nencini, "Hypercomplex
%                            quality assessment of multi/hyper-spectral
%                            images," IEEE Geoscience and Remote Sensing
%                            Letters, vol. 6, pp. 662–665, Oct. 2009.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Q2n_index,Q2n_map] = q2nf(I3,X3,bs)

[N1,N2,N3]=size(I3);

stepx=ceil(N1/bs);
stepy=ceil(N2/bs);

if stepy<=0
    stepy=1;
    stepx=1;
end;

est1=(stepx-1)*bs+bs-N1;
est2=(stepy-1)*bs+bs-N2;

if sum([(est1~=0),(est2~=0)])>0
  refref=[];
  fusfus=[];
  
  for ii=1:N3
      a1=squeeze(I3(:,:,1));
    
      ia1=zeros(N1+est1,N2+est2);
      ia1(1:N1,1:N2)=a1;
      ia1(:,N2+1:N2+est2)=ia1(:,N2:-1:N2-est2+1);
      ia1(N1+1:N1+est1,:)=ia1(N1:-1:N1-est1+1,:);
      refref=cat(3,refref,ia1);
      
      if ii<N3
          I3=I3(:,:,2:end);
      end
  end

  I3=refref;
  clear refref
  
  for ii=1:N3
      a2=squeeze(X3(:,:,1));
      
      ia2=zeros(N1+est1,N2+est2);
      ia2(1:N1,1:N2)=a2;
      ia2(:,N2+1:N2+est2)=ia2(:,N2:-1:N2-est2+1);
      ia2(N1+1:N1+est1,:)=ia2(N1:-1:N1-est1+1,:);
      fusfus=cat(3,fusfus,ia2);
      
      if ii<N3
          X3=X3(:,:,2:end);
      end
  end
  
  X3=fusfus;
  clear fusfus a1 a2 ia1 ia2

end

%X3=uint16(X3);
%I3=uint16(I3);

[N1,N2,N3]=size(I3);

if ceil(log2(N3))-log2(N3)~=0
    Ndif=2^(ceil(log2(N3)))-N3;
    dif=zeros(N1,N2,Ndif);
    %dif=uint16(dif);
    I3=cat(3,I3,dif);
    X3=cat(3,X3,dif);
end
[~,~,N3]=size(I3);

valori=zeros(stepx,stepy,N3);

for jj=1:stepx
    for ii=1:stepy
        o=onions_quality(I3((jj-1)*bs+1:(jj-1)*bs+bs,(ii-1)*bs+1:(ii-1)*bs+bs,:),...
                         X3((jj-1)*bs+1:(jj-1)*bs+bs,(ii-1)*bs+1:(ii-1)*bs+bs,:),bs);
        valori(jj,ii,:)=o;    
    end
end

Q2n_map=sqrt(sum((valori.^2),3));

Q2n_index=mean2(Q2n_map);

end