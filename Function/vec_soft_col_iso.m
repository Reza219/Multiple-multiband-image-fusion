function [Y1,Y2] = vec_soft_col_iso(X1,X2,t)
%  computes the isotropic vector soft columnwise

n = sqrt(sum(X1.^2)+sum(X2.^2));
a = max(n-t,0);
A = repmat((a./n),size(X1,1),1);
Y1 = A.*X1;
Y2 = A.*X2;

end