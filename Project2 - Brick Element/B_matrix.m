function [B] = B_matrix( dNdx, dNdy, dNdz )
% Preparing B matrix
Bi=zeros(6,3);
for i=1:8 
    % Find Bi at each node
    Bi(i)=[dNdx(i) 0       0;
           0       dNdy(i) 0;
           0       0       dNdz(i);
           dNdy(i) dNdx(i) 0;
           0       dNdz(i) dNdy(i);
           dNdz(i) 0       dNdx(i)];
end
% B matrix
B=[Bi(1) Bi(2) Bi(3) Bi(4) Bi(5) Bi(6) Bi(7) Bi(8)]; % Assembles 6x24 B matrix
end