function [Jac,dNdx,dNdy,dNdz] = Jacobian(Xi_N,Eta_N,Zeta_N,Xi,Eta,Zeta,x,y,z)
% This finds Jacobian at current value of Xi,Eta and Zeta
% Also, it finds the derivatives of shape functions w.r.t x,y and z

% Making zero matrices for the derivative of physical and natural coordinates
dNdXi=zeros(8,1); dNdEta=zeros(8,1); dNdZeta=zeros(8,1); 
dNdx=zeros(8,1); dNdy=zeros(8,1); dNdz=zeros(8,1);
    
for i=1:8
% Finding the derivative of shape function w.r.t physical coordinates
    dNdXi(i)=(1/8)*Xi_N(i)*(1+Eta*Eta_N(i))*(1+Zeta*Zeta_N(i));
    dNdEta(i)=(1/8)*Eta_N(i)*(1+Xi*Xi_N(i))*(1+Zeta*Zeta_N(i));
    dNdZeta(i)=(1/8)*Zeta_N(i)*(1+Xi*Xi_N(i))*(1+Eta*Eta_N(i));
end

% Finding a Jacobian
Jac=[dNdXi dNdEta dNdZeta]'*[x' y' z']; % Jacobian [3x3]

for i=1:8
% The derivative of N w.r.t physical coordinates
    dN=Jac\[dNdXi(i) dNdEta(i) dNdZeta(i)]'; 
    dNdx(i)=dN(1); % the derivative of shape function w.r.t x
    dNdy(i)=dN(2); % the derivative of shape function w.r.t y
    dNdz(i)=dN(3); % the derivative of shape function w.r.t z
end
end 