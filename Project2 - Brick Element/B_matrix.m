function [B] = B_matrix(dNdx,dNdy,dNdz)
% Preparing B matrix [6x24]
B=zeros(6,24); 

for i=1:8 
% Finding Bi at each node, there are 8 of it
    Bi=[dNdx(i)      0        0;
           0       dNdy(i)    0;
           0         0      dNdz(i);
        dNdy(i)    dNdx(i)    0;
           0       dNdz(i)  dNdy(i);
        dNdz(i)      0      dNdx(i)];
% B matrix 
B(:,(i*3-2):(i*3))=Bi(:,:); % This step is for assembling a 6x24 B matrix
end
end