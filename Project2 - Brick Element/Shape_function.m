function [ N ] = Shape_function(Xi_I,Eta_I,Zeta_I,Xi,Eta,Zeta)
% Shape functions for brick element. 

N=zeros(3,24);
for i=1:8
    Ni=(1/8)*(1+Xi*Xi_I(i))*(1+Eta*Eta_I(i))*(1+Zeta*Zeta_I(i)); % it finds shape function for each node
    Nd=[Ni 0 0;
        0 Ni 0;
        0 0 Ni]; % Creates diagonal matrix to create N matrix later
    
    % Shape function, N
    N(:,(i*3-2):(i*3))=Nd(:,:); % This assembles 3x24 N matrix 
end
end