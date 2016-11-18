function [ Emat ] = E_matrix(E,G)
% Material Matrix

nu=(E/(2*G))-1; % Poisson's Ratio
nu1=1-nu; 
nu2=1+nu;
nu3=1-2*nu;
Emat=zeros(6,6); % Creates a matrix of all zeros
Emat(1:3,1:3)=nu; % Fills the top left 3x3 block with nu
Emat(1,1)=nu1; % Overwrites the diagonal
Emat(2,2)=nu1;
Emat(3,3)=nu1;
Emat(4,4)=nu3/2; % Writes the bottom right diagonal
Emat(5,5)=nu3/2; 
Emat(6,6)=nu3/2;
Emat=Emat*E/(nu2*nu3); % Multiplies by the outside constants

end
