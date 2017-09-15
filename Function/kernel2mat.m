function B=kernel2mat(K,nr,nc)

% flip the kernel
K=rot90(K,2);
mc=round((nc+1)/2);
mr=round((nr+1)/2);

% the size of the kernel
[kh,kv]=size(K);
lx=(kh-1)/2;
ly=(kv-1)/2;
B=zeros(nr,nc);

% range of the pixels
B(mr-lx:mr+lx,mc-ly:mc+ly)=K;
B=circshift(B,[-mr+1,-mc+1]);