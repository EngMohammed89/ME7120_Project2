function [ J,dNdx,dNdy,dNdz ] = J_brick( Epsilon, Eta, Zeta, x, y, z )
% J matrix. 
% Epsilon, Eta, Zeta are determined by the location in the brick element
EpsilonI=[-1 -1 -1 -1 1 1 1 1]; 
EtaI=[-1 -1 1 1 -1 -1 1 1];
ZetaI=[1 -1 -1 1 1 -1 -1 1];

for i=1:8 % iterates over each node
    dNdEpsilon(i)=(1/8)*EpsilonI(i)*(1+Eta*EtaI(i))*(1+Zeta*ZetaI(i));
    dNdEta(i)=(1/8)*EtaI(i)*(1+Epsilon*EpsilonI(i))*(1+Zeta*ZetaI(i));
    dNdZeta(i)=(1/8)*ZetaI(i)*(1+Epsilon*EpsilonI(i))*(1+Eta*EtaI(i));
end

J=[dNdEpsilon; dNdEta; dNdZeta]*[x' y' z']; 
Jp=J^-1; 

for i=1:8
    Direction=Jp*[dNdEpsilon(i); dNdEta(i); dNdZeta(i)]; % finds the cartesian derivatives
    dNdx(i)=Direction(1); % stores the current nodal x derivative
    dNdy(i)=Direction(2); % stores the current nodal y derivative
    dNdz(i)=Direction(3); % stores the current nodal z derivative
end
