function [ B ] = B_matrix(dNdx,dNdy,dNdz)
% Preparing B matrix

B=zeros(6,24); 
for i=1:8 
    % Find Bi at each node
    Bi=[dNdx(i)      0        0;
           0       dNdy(i)    0;
           0         0      dNdz(i);
        dNdy(i)    dNdx(i)    0;
           0       dNdz(i)  dNdy(i);
        dNdz(i)      0      dNdx(i)];
% B matrix 
B(:,(i*3-2):(i*3))=Bi(:,:); % This assembles 6x24 B matrix
end
end
