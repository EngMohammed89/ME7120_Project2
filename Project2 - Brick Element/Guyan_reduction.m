function [ Kr,Mr ] = Guyan_reduction( Ke,Me )
% Guyan reduction is used to reduce the size of stiffness and mass matrix 
% form 33x33 to 24x24 

K11=Ke(1:24,1:24);
K21=Ke(25:33,1:24);
K22=Ke(25:33,25:33);
I=eye(size(K11,1));
Q =[I;
   -K22\K21];
% Reduced stiffness matrix
Kr=Q'*Ke*Q;
% Reduced mass matrix
Mr=Me;
end
