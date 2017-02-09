function [Y]=PPlus_s(X,n_dr,n_dc)
% when input is 
%      X=[X11 X12 X13 X14;
%         X21 X22 X23 X24;
%         X31 X32 X33 X34;
%         X41 X42 X43 X44]
% the output is 
%      X=[\sum_{i,j}X_{i,j} X12 X13 X14;
%         X21 X22 X23 X24;
%         X31 X32 X33 X34;
%         X41 X42 X43 X44]

[nr,nc,nb]=size(X);

%% Sum according to the column
Temp=reshape(X,[nr*n_dc nc/n_dc nb]);
Temp(:,1,:)=sum(Temp,2);

%% Sum according to the row
Temp1=reshape(permute(reshape(Temp(:,1,:),[nr n_dc nb]),[2 1 3]),[n_dc*n_dr nr/n_dr nb]);
Y=permute(reshape(sum(Temp1,2),[n_dc n_dr nb]),[2 1 3]);

end
