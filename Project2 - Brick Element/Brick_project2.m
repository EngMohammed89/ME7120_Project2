function out=Brick_project2(mode,b,c,d,e)
  % Brick element does as listed below. 
% Its properties (bprops) are in the order
% bprops=[E G rho]
% See wfem.m for more explanation.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Variables (global):
% -------------------
% K       :    Global stiffness matrix
% Ks      :    Global stiffness buckling matrix
% M       :    Global mass matrix
% nodes   :    [x y z] nodal locations

global ismatnewer
global K
global Ks
global M
global nodes % Node locations
global elprops
global element
global points
global Fepsn % Initial strain "forces". 
global lines
global restart
global reload
global curlineno
global DoverL
global surfs
%
% Variables (local):
% ------------------
% bnodes  :    node/point numbers for actual brick nodes 1-8
% k       :    stiffness matrix in local coordiates
% kg      :    stiffness matrix rotated into global coordinates
% m       :    mass matrix in local coordiates
% mg      :    mass matrix rotated into global coordinates
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

out=0;
if strcmp(mode,'numofnodes')
    % This allows a code to find out how many nodes this element has
    out=8;
end
if strcmp(mode,'generate')
  elnum=c;%When this mode is called, the element number is the 3rd
          %argument. 
  
          %The second argument (b) is the element
          %definition. For this element b is
          %node1 node2 node3 node4 node5 node6 node7 node8 materialnumber
  
          %There have to be 9 elements for this element's
          %definition (above)
  if length(b)==9
      element(elnum).nodes=b(1:8);
      element(elnum).properties=b(9);
  else 
	  b
      %There have to be 9 numbers on a line defining the element. 
      warndlg(['Element ' num2str(elnum) ' on line ' ...
               num2str(element(elnum).lineno) ' entered incorrectly.'], ...
              ['Malformed Element'],'modal')
      return
  end
 
end

if strcmp(mode,'make')||strcmp(mode,'istrainforces')
  elnum=b;% When this mode is called, the element number is given
          % as the second input.
bnodes=[element(elnum).nodes];    
bprops=elprops(element(elnum).properties).a;% element(elnum).properties 
                                              % stores the
                                              % properties number
                                              % of the current
                                              % elnum. elprops
                                              % contains this
                                              % data. This is
                                              % precisely the
                                              % material properties
                                              % line in an
                                              % array. You can pull
                                              % out any value you
                                              % need for your use. 
  
% 
  if length(bprops)==3
      Em=bprops(1);  % Modulus of elasticity
      G=bprops(2);
      rho=bprops(3);
  else
      warndlg(['The number of material properties set for ' ...
               'this element (' num2str(length(bprops)) ') isn''t ' ...
               'appropriate for a beam3 element. '      ...
               'Please refer to the manual.'],...
              'Bad element property definition.','modal');
  end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Brick properties (bprops) are in the order
% bprops=[E G rho]


if strcmp(mode,'make')

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  % Define beam node locations for easy later referencing
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  x1=nodes(bnodes(1),1);
  y1=nodes(bnodes(1),2);
  z1=nodes(bnodes(1),3);
  x2=nodes(bnodes(2),1);
  y2=nodes(bnodes(2),2);
  z2=nodes(bnodes(2),3);
  x3=nodes(bnodes(3),1);
  y3=nodes(bnodes(3),2);
  z3=nodes(bnodes(3),3);
  x4=nodes(bnodes(4),1);
  y4=nodes(bnodes(4),2);
  z4=nodes(bnodes(4),3);
  x5=nodes(bnodes(5),1);
  y5=nodes(bnodes(5),2);
  z5=nodes(bnodes(5),3);
  x6=nodes(bnodes(6),1);
  y6=nodes(bnodes(6),2);
  z6=nodes(bnodes(6),3);
  x7=nodes(bnodes(7),1);
  y7=nodes(bnodes(7),2);
  z7=nodes(bnodes(7),3);
  x8=nodes(bnodes(8),1);
  y8=nodes(bnodes(8),2);
  z8=nodes(bnodes(8),3);
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 %% x, y and z are nodal locations. They are used to find Jacobian later.
  x=[x1 x2 x3 x4 x5 x6 x7 x8];
  y=[y1 y2 y3 y4 y5 y6 y7 y8];
  z=[z1 z2 z3 z4 z5 z6 z7 z8]; 
  
 %% The nodal values of the natural coordinates
  Xi_N=[-1 -1 -1 -1 1 1 1 1]; 
  Eta_N=[-1 -1 1 1 -1 -1 1 1];
  Zeta_N=[1 -1 -1 1 1 -1 -1 1];
  
 %% Elasticity Matrix  
v = (Em/(2*G))-1;
E = (Em/((1+v)*(1-2*v)))* [(1-v)  v    v      0        0        0;
                             v  (1-v)  v      0        0        0;
                             v    v  (1-v)    0        0        0;
                             0    0    0  (1-2*v)/2    0        0;
                             0    0    0      0    (1-2*v)/2    0;
                             0    0    0      0        0   (1-2*v)/2];

  %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  % Define gauss points and loop to integrate for Ke and Me
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  Me=zeros(24,24); 
  Ke=zeros(24,24); 
  
  % Derivation of the stiffness matrix
  numbrickgauss=5; % Number of gauss points for the stiffness matrix
  [bgpts,bgpw]=gauss([numbrickgauss,numbrickgauss,numbrickgauss]); % The gauss points and weights
  for i=1:size(bgpts,1) 
      [Jac,dNdx,dNdy,dNdz]=Jacobian(Xi_N,Eta_N,Zeta_N,bgpts(i,1),bgpts(i,2),bgpts(i,3),x,y,z); % Jacobian for the current Gauss point
      B=B_matrix(dNdx,dNdy,dNdz); % B matrix for the current Gauss point    
      Ke=Ke+bgpw(i)*B'*E*B*det(Jac); % Gauss integration part to find the stiffness matrix
  end
  
  % Derivation of the mass matrix
  numbrickgauss=numbrickgauss+2; % number of gauss points for the mass matrix
  [bgpts,bgpw]=gauss([numbrickgauss,numbrickgauss,numbrickgauss]); % the gauss points and weights
  for i=1:size(bgpts,1) 
      Jac=Jacobian(Xi_N,Eta_N,Zeta_N,bgpts(i,1),bgpts(i,2),bgpts(i,3),x,y,z); % Jacobian for the current Gauss point
      N=Shape_function(Xi_N,Eta_N,Zeta_N,bgpts(i,1),bgpts(i,2),bgpts(i,3)); % Shape function N for the current Gauss point
      Me=Me+bgpw(i)*N'*rho*N*det(Jac); % Gauss integration part to find the element mass matrix
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  % Assembling matrices into global matrices
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  bn1=bnodes(1);bn2=bnodes(2);bn3=bnodes(3);bn4=bnodes(4);
  bn5=bnodes(5);bn6=bnodes(6);bn7=bnodes(7);bn8=bnodes(8);
  
  indices=[bn1*6+(-5:-3) bn2*6+(-5:-3) bn3*6+(-5:-3) bn4*6+(-5:-3)...
           bn5*6+(-5:-3) bn6*6+(-5:-3) bn7*6+(-5:-3) bn8*6+(-5:-3)] ;


  K(indices,indices)=K(indices,indices)+Ke;
  M(indices,indices)=M(indices,indices)+Me;

  % At this point we also know how to draw the element (what lines
  % and surfaces exist). For the beam3 element, 2 lines are
  % appropriate. Just add the pair of node numbers to the lines
  % array and that line will always be drawn.
  numlines=size(lines,1);
  lines(numlines+1,:)=[bn1 bn2];
  lines(numlines+2,:)=[bn2 bn3];
  lines(numlines+3,:)=[bn3 bn4];
  lines(numlines+4,:)=[bn4 bn1];
  lines(numlines+5,:)=[bn5 bn6];
  lines(numlines+6,:)=[bn6 bn7];
  lines(numlines+7,:)=[bn7 bn8];
  lines(numlines+8,:)=[bn8 bn5];
  lines(numlines+9,:)=[bn1 bn5];
  lines(numlines+10,:)=[bn2 bn6];
  lines(numlines+11,:)=[bn3 bn7];
  lines(numlines+12,:)=[bn4 bn8];
  
  
  panelcolor=[1 0 1];% This picks a color. You can change the
                     % numbes between 0 and 1. 
  %Don't like this color? Use colorui to pick another one. Another
  %option is that if we can't see the elements separately we can
  %chunk up x*y*z, divide by x*y*x of element, see if we get
  %integer powers or not to define colors that vary by panel. 
  
  panelcolor2=[0 1 1];panelcolor3=[1 0 0];panelcolor4=[1 1 1]; panelcolor5=[0 0 0];panelcolor6=[0 1 0];
  surfs=[surfs;bn1 bn2 bn3 bn4 panelcolor];
  surfs=[surfs;bn2 bn6 bn7 bn3 panelcolor2];
  surfs=[surfs;bn6 bn5 bn8 bn7 panelcolor3];
  surfs=[surfs;bn5 bn1 bn4 bn8 panelcolor4];
  surfs=[surfs;bn1 bn5 bn6 bn2 panelcolor5];
  surfs=[surfs;bn4 bn8 bn7 bn3 panelcolor6];

  
  %Each surface can have a different color if you like. Just change
  %the last three numbers on the row corresponding to that
  %surface. 

elseif strcmp(mode,'istrainforces')
elseif strcmp(mode,'draw')
elseif strcmp(mode,'buckle')
end
