function [ N ] = Shape_function(Xi_N,Eta_N,Zeta_N,Xi,Eta,Zeta)
% Shape functions, N [3x24] 

N=zeros(3,24);
for i=1:8
    Ni=(1/8)*(1+Xi*Xi_N(i))*(1+Eta*Eta_N(i))*(1+Zeta*Zeta_N(i)); % it finds shape function for each node
% creating a diagonal matrix to assemble N matrix later    
    Nd=[Ni 0 0;
        0 Ni 0;
        0 0 Ni];    
% Shape function, N
N(:,(i*3-2):(i*3))=Nd(:,:); % assembling a [3x24] N matrix 
end
end