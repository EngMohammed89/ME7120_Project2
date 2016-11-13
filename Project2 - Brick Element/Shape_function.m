function [ N ] = Shape_function(Xi_I,Eta_I,Zeta_I,Xi,Eta,Zeta)
% Shape functions for brick element. 
Nd=zeros(3,3)
for i=1:8 % iterates over each node
    Ni=(1/8)*(1+Xi*Xi_I(i))*(1+Eta*Eta_I(i))*(1+Zeta*Zeta_I(i)); % it finds shape function for each node
    Nd(i)=[Ni 0 0;
           0 Ni 0;
           0 0 Ni]; % Creates diagonal matrix to create N matrix later
end
% Shape function N
N=[Nd(1) Nd(2) Nd(3) Nd(4) Nd(5) Nd(6) Nd(7) Nd(8)]; % Assembles 6x24 B matrix
end