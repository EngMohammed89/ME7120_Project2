function [ Jac,dNdx,dNdy,dNdz ] = Jacobian(Xi_I,Eta_I,Zeta_I,Xi,Eta,Zeta,x,y,z)
% This finds Jacobian at current value of Xi,Eta and Zeta
% Also, it finds the derivatives of shape functions w.r.t x,y and z


% Makes zero matrices for the derivative of physical and natural coordinates
dNdXi=zeros(8,1); dNdEta=zeros(8,1); dNdZeta=zeros(8,1); 
dNdx=zeros(8,1); dNdy=zeros(8,1); dNdz=zeros(8,1);
    
for i=1:8
    dNdXi(i)=(1/8)*Xi_I(i)*(1+Eta*Eta_I(i))*(1+Zeta*Zeta_I(i));
    dNdEta(i)=(1/8)*Eta_I(i)*(1+Xi*Xi_I(i))*(1+Zeta*Zeta_I(i));
    dNdZeta(i)=(1/8)*Zeta_I(i)*(1+Xi*Xi_I(i))*(1+Eta*Eta_I(i));
end

Jac=[dNdXi dNdEta dNdZeta]'*[x' y' z']; % Jacobian [3x3]

for i=1:8
R=Jac\[dNdXi(i) dNdEta(i) dNdZeta(i)]'; % finds the cartesian derivatives
    dNdx(i)=R(1); % the derivative of shape function w.r.t x
    dNdy(i)=R(2); % the derivative of shape function w.r.t y
    dNdz(i)=R(3); % the derivative of shape function w.r.t z
end
end 