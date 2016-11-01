clear all
%% Pyramid Brick Element

% locations of 8 nodes of the brick in XYZ co-ordinate system
X1 = 0; X2 = 0; X3 = 32; X4 = 32; X5 = 74; X6 = 74; X7 = 42; X8 = 42; 

Y1 = 0; Y2 = 0; Y3 = 190; Y4 = 190; Y5 = 0; Y6 = 0; Y7 = 190; Y8 = 190;

Z1 = 74; Z2 = 0; Z3 = 32; Z4 = 42; Z5 = 74; Z6 = 0; Z7 = 32; Z8 = 42;

% Natural co-ordinates for Jacobian at middle
Xi = -1
Eta = -1
Zeta = 1

% Shape Funtions for Brick element (put co-ordinate values for Xi, Eta and
% Zeta in (i) in Ni = 1/8 (1+Xi*Xi(i))*(1+Eta*Eta(i))*(1+Zeta*Zeta(i))
N1 = 1/8*(1-Xi)*(1-Eta)*(1+Zeta);
N2 = 1/8*(1-Xi)*(1-Eta)*(1-Zeta);
N3 = 1/8*(1-Xi)*(1+Eta)*(1-Zeta);
N4 = 1/8*(1-Xi)*(1+Eta)*(1+Zeta);
N5 = 1/8*(1+Xi)*(1-Eta)*(1+Zeta);
N6 = 1/8*(1+Xi)*(1-Eta)*(1-Zeta);
N7 = 1/8*(1+Xi)*(1+Eta)*(1-Zeta);
N8 = 1/8*(1+Xi)*(1+Eta)*(1+Zeta);

% Assmbly of [3x3] Jacobian matrix 
J11 = (X1*(Eta - 1)*(Zeta + 1))/8 - (X2*(Eta - 1)*(Zeta - 1))/8 + (X3*(Eta + 1)*(Zeta - 1))/8 - (X4*(Eta + 1)*(Zeta + 1))/8 - (X5*(Eta - 1)*(Zeta + 1))/8 + (X6*(Eta - 1)*(Zeta - 1))/8 - (X7*(Eta + 1)*(Zeta - 1))/8 + (X8*(Eta + 1)*(Zeta + 1))/8;
J12 = (Y1*(Eta - 1)*(Zeta + 1))/8 - (Y2*(Eta - 1)*(Zeta - 1))/8 + (Y3*(Eta + 1)*(Zeta - 1))/8 - (Y4*(Eta + 1)*(Zeta + 1))/8 - (Y5*(Eta - 1)*(Zeta + 1))/8 + (Y6*(Eta - 1)*(Zeta - 1))/8 - (Y7*(Eta + 1)*(Zeta - 1))/8 + (Y8*(Eta + 1)*(Zeta + 1))/8;
J13 = (Z1*(Eta - 1)*(Zeta + 1))/8 - (Z2*(Eta - 1)*(Zeta - 1))/8 + (Z3*(Eta + 1)*(Zeta - 1))/8 - (Z4*(Eta + 1)*(Zeta + 1))/8 - (Z5*(Eta - 1)*(Zeta + 1))/8 + (Z6*(Eta - 1)*(Zeta - 1))/8 - (Z7*(Eta + 1)*(Zeta - 1))/8 + (Z8*(Eta + 1)*(Zeta + 1))/8;
J21 = X1*(Xi/8 - 1/8)*(Zeta + 1) - X2*(Xi/8 - 1/8)*(Zeta - 1) + X3*(Xi/8 - 1/8)*(Zeta - 1) - X4*(Xi/8 - 1/8)*(Zeta + 1) - X5*(Xi/8 + 1/8)*(Zeta + 1) + X6*(Xi/8 + 1/8)*(Zeta - 1) - X7*(Xi/8 + 1/8)*(Zeta - 1) + X8*(Xi/8 + 1/8)*(Zeta + 1);
J22 = Y1*(Xi/8 - 1/8)*(Zeta + 1) - Y2*(Xi/8 - 1/8)*(Zeta - 1) + Y3*(Xi/8 - 1/8)*(Zeta - 1) - Y4*(Xi/8 - 1/8)*(Zeta + 1) - Y5*(Xi/8 + 1/8)*(Zeta + 1) + Y6*(Xi/8 + 1/8)*(Zeta - 1) - Y7*(Xi/8 + 1/8)*(Zeta - 1) + Y8*(Xi/8 + 1/8)*(Zeta + 1);
J23 = Z1*(Xi/8 - 1/8)*(Zeta + 1) - Z2*(Xi/8 - 1/8)*(Zeta - 1) + Z3*(Xi/8 - 1/8)*(Zeta - 1) - Z4*(Xi/8 - 1/8)*(Zeta + 1) - Z5*(Xi/8 + 1/8)*(Zeta + 1) + Z6*(Xi/8 + 1/8)*(Zeta - 1) - Z7*(Xi/8 + 1/8)*(Zeta - 1) + Z8*(Xi/8 + 1/8)*(Zeta + 1);
J31 = X1*(Xi/8 - 1/8)*(Eta - 1) - X2*(Xi/8 - 1/8)*(Eta - 1) + X3*(Xi/8 - 1/8)*(Eta + 1) - X4*(Xi/8 - 1/8)*(Eta + 1) - X5*(Xi/8 + 1/8)*(Eta - 1) + X6*(Xi/8 + 1/8)*(Eta - 1) - X7*(Xi/8 + 1/8)*(Eta + 1) + X8*(Xi/8 + 1/8)*(Eta + 1);
J32 = Y1*(Xi/8 - 1/8)*(Eta - 1) - Y2*(Xi/8 - 1/8)*(Eta - 1) + Y3*(Xi/8 - 1/8)*(Eta + 1) - Y4*(Xi/8 - 1/8)*(Eta + 1) - Y5*(Xi/8 + 1/8)*(Eta - 1) + Y6*(Xi/8 + 1/8)*(Eta - 1) - Y7*(Xi/8 + 1/8)*(Eta + 1) + Y8*(Xi/8 + 1/8)*(Eta + 1);
J33 = Z1*(Xi/8 - 1/8)*(Eta - 1) - Z2*(Xi/8 - 1/8)*(Eta - 1) + Z3*(Xi/8 - 1/8)*(Eta + 1) - Z4*(Xi/8 - 1/8)*(Eta + 1) - Z5*(Xi/8 + 1/8)*(Eta - 1) + Z6*(Xi/8 + 1/8)*(Eta - 1) - Z7*(Xi/8 + 1/8)*(Eta + 1) + Z8*(Xi/8 + 1/8)*(Eta + 1);


% Jacobian at middle
J = [J11 J12 J13;
     J21 J22 J23;
     J31 J32 J33];
 
% Partial derivatives shape functions for [B] w.r.t. X, Y and Z
C1 = inv(J)*  [((Eta - 1)*(Zeta + 1))/8; 
                 (Xi/8 - 1/8)*(Zeta + 1); 
                  (Xi/8 - 1/8)*(Eta - 1)];
       
B1 = [C1(1)     0         0;
      0       C1(2)       0;
      0         0     C1(3);
      C1(2)   C1(2)       0;
      0       C1(3)   C1(2);    
      C1(3)     0     C1(1)];
  
C2 = inv(J)*  [-((Eta - 1)*(Zeta - 1))/8; 
                 -(Xi/8 - 1/8)*(Zeta - 1); 
                  -(Xi/8 - 1/8)*(Eta - 1)];
       
B2 = [C2(1)     0         0;
      0       C2(2)       0;
      0         0     C2(3);
      C2(2)   C2(2)       0;
      0       C2(3)   C2(2);    
      C2(3)     0     C2(1)];
  
C3 = inv(J)*  [((Eta + 1)*(Zeta - 1))/8; 
                 (Xi/8 - 1/8)*(Zeta - 1); 
                  (Xi/8 - 1/8)*(Eta + 1)];
       
B3 = [C3(1)     0         0;
      0       C3(2)       0;
      0         0     C3(3);
      C3(2)   C3(2)       0;
      0       C3(3)   C3(2);    
      C3(3)     0     C3(1)]; 

C4 = inv(J)*  [-((Eta + 1)*(Zeta + 1))/8; 
                 -(Xi/8 - 1/8)*(Zeta + 1); 
                 -(Xi/8 - 1/8)*(Eta + 1)];
       
B4 = [C4(1)     0         0;
      0       C4(2)       0;
      0         0     C4(3);
      C4(2)   C4(2)       0;
      0       C4(3)   C4(2);    
      C4(3)     0     C4(1)];
 
C5 = inv(J)*  [-((Eta - 1)*(Zeta + 1))/8; 
                 -(Xi/8 + 1/8)*(Zeta + 1); 
                 -(Xi/8 + 1/8)*(Eta - 1)];
       
B5 = [C5(1)     0         0;
      0       C5(2)       0;
      0         0     C5(3);
      C5(2)   C5(2)       0;
      0       C5(3)   C5(2);    
      C5(3)     0     C5(1)]; 

C6 = inv(J)*  [((Eta - 1)*(Zeta - 1))/8; 
                 (Xi/8 + 1/8)*(Zeta - 1); 
                 (Xi/8 + 1/8)*(Eta - 1)];
       
B6 = [C6(1)     0         0;
      0       C6(2)       0;
      0         0     C6(3);
      C6(2)   C6(2)       0;
      0       C6(3)   C6(2);    
      C6(3)     0     C6(1)];
  
C7 = inv(J)*  [-((Eta + 1)*(Zeta - 1))/8; 
                 -(Xi/8 + 1/8)*(Zeta - 1); 
                 -(Xi/8 + 1/8)*(Eta + 1)];
       
B7 = [C7(1)     0         0;
      0       C7(2)       0;
      0         0     C7(3);
      C7(2)   C7(2)       0;
      0       C7(3)   C7(2);    
      C7(3)     0     C7(1)];   

C8 = inv(J)*  [((Eta + 1)*(Zeta + 1))/8; 
                 (Xi/8 + 1/8)*(Zeta + 1); 
                 (Xi/8 + 1/8)*(Eta + 1)];
       
B8 = [C8(1)     0         0;
      0       C8(2)       0;
      0         0     C8(3);
      C8(2)   C8(2)       0;
      0       C8(3)   C8(2);    
      C8(3)     0     C8(1)];   

B = [B1 B2 B3 B4 B5 B6 B7 B8];
  
% Determining Elasticity matrix [E] 
% Steel Properties
modulus = 203.4e9;
v = 0.3;
  
E = (modulus/((1+v)*(1-2*v)))* [(1-v)  v    v      0        0        0;
                                  v  (1-v)  v      0        0        0;
                                  v    v  (1-v)    0        0        0;
                                  0    0    0  (1-2*v)/2    0        0;
                                  0    0    0      0    (1-2*v)/2    0;
                                  0    0    0      0        0   (1-2*v)/2];
                              

Ke = B' * E * det(J) * 8                          
                              
                              
                              
                              
                              
      