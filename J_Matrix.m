function [ J,dNdx,dNdy,dNdz ] = J_brick(Xi, Eta, Zeta, x, y, z)
% J Matrix based on nodal locations in brick of xi, eta, and zeta

XiI=   [-1 -1 -1 -1  1  1  1  1];
EtaI=  [-1 -1  1  1 -1 -1  1  1];
ZetaI= [-1  1  1 -1 -1  1  1 -1];

for i=1:8 
    dNdEpsilon(i)=(1/8)*XiI(i)*(1+Eta*EtaI(i))*(1+Zeta*ZetaI(i));
    dNdEta(i)=(1/8)*EtaI(i)*(1+Xi*XiI(i))*(1+Zeta*ZetaI(i));
    dNdZeta(i)=(1/8)*Zeta(i)*(1+Xi*XiI(i))*(1+Eta*EtaI(i));
end

J=[dNdXi; dNdEta; dNdZeta]*[x' y' z']; % assembles derivative and nodal locations matrices and multiplies
Jp=J^-1;

for i=1:8
    Resultant=Jp*[dNdXi(i); dNdEta(i); dNdZeta(i)]; %
    dNdx(i)=Resultant(1); 
    dNdy(i)=Resultant(2); 
    dNdz(i)=Resultant(3); 
end
