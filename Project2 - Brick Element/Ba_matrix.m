function [ Ba ] = Ba_matrix(J0,Xi,eta,zeta)
% Preparing Ba matrix, Ba=alpha*beta*gamma
% finding alpha
alpha=zeros(6,9);
alpha(1,1)=1;alpha(2,5)=1;alpha(3,9)=1;
alpha(4,2)=1;alpha(4,4)=1;alpha(5,6)=1;
alpha(5,8)=1;alpha(6,3)=1;alpha(6,7)=1;

% finding beta
beta=zeros(9,9);I=eye(3,3)
J_inv=J0\I;
beta(1:3,1:3)=J_inv;beta(4:6,4:6)=J_inv;
beta(7:9,7:9)=J_inv;

% finding gamma
gamma=zeros(9,9);
gamma(1,1)=-2*Xi; gamma(4,2)=-2*Xi; gamma(7,3)=-2*Xi;
gamma(2,4)=-2*eta; gamma(5,5)=-2*eta; gamma(8,6)=-2*eta;
gamma(3,7)=-2*zeta; gamma(6,8)=-2*zeta; gamma(9,9)=-2*zeta;

Ba=alpha*beta*gamma
end

