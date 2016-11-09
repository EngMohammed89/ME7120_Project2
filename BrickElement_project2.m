function out=BrickElement_project2(mode,b,c,d,e)
  
% Brick element does as listed below. 
% Its properties (bprops) are in the order
% bprops=[E G rho v]
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
% bnodes  :    node/point numbers for actual beam nodes 1-2 and point
% k       :    stiffness matrix in local coordiates
% kg      :    stiffness matrix rotated into global coordinates
% m       :    mass matrix in local coordiates
% mg      :    mass matrix rotated into global coordinates

out=0;
if strcmp(mode,'numofnodes')
    % This allows a code to find out how many nodes this element has
    out=8;
end
if strcmp(mode,'generate')
  elnum=c;%When this mode is called, 
          %The argument (b) is the element
          %definition. For this element b is
          %node1-node8 point(for rotation) and material#
          %There have to be 10 elements for this element's
          %definition (above)
          
  if length(b)==10
      element(elnum).nodes=b(1:8);
      element(elnum).properties=b(10);
      element(elnum).point=b(9);
  else 
	  b
      %There have to be ten numbers on a line defining the
      %element. 
      warndlg(['Element ' num2str(elnum) ' on line ' ...
               num2str(element(elnum).lineno) ' entered incorrectly.'], ...
              ['Malformed Element'],'modal')
      return
  end
 
end

% Here we figure out what the beam properties mean. If you need
% them in a mode, that mode should be in the if on the next line.
if strcmp(mode,'make')||strcmp(mode,'istrainforces')
  elnum=b;% When this mode is called, the element number is given
          % as the second input.
  bnodes=[element(elnum).nodes element(elnum).point];% The point is
                                                     % referred to
                                                     % as node 4
                                                     % below,
                                                     % although it
                                                     % actually
                                                     % calls the
                                                     % array points
                                                     % to get its
                                                     % location. Its
                                                     % not really a
                                                     % node, but
                                                     % just a point
                                                     % that helps
                                                     % define
                                                     % orientation. Your
                                                     % element may
                                                     % not need
                                                     % such a
                                                     % reference point.
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
  if length(bprops)==4
      Em=bprops(1); % Modulus of elasticity
      G=bprops(2);
      rho=bprops(3);
      v=bprops(4);
%       A1=bprops(4);
%       A2=bprops(5);
%       J1=bprops(6);
%       J2=bprops(7);
%       Izz1=bprops(8);
%       Izz2=bprops(9);
%       Iyy1=bprops(10);
%       Iyy2=bprops(11);
  else
      warndlg(['The number of material properties set for ' ...
               'this element (' num2str(length(bprops)) ') isn''t ' ...
               'appropriate for a Brick Element. '      ...
               'Please refer to the manual.'],...
              'Bad element property definition.','modal');
  end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Brick element properties (bprops) are in the order
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
  x9=points(bnodes(9),1);
  y9=points(bnodes(9),2);
  z9=points(bnodes(9),3);
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  % Shape functions for higher order beam. 
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Shape functions in matrix polynomial form (polyval style) for bending
  % Should be 4 shape functions, each a 1 x 4 matrix
  bn1 =   [0.25 0.00 -0.75 0.50];
  bn1d =  [0.75 0.00 -0.75];
  bn1dd = [1.50 0.00];
  bn2 =   [0.25 -0.25 -0.25 0.25];
  bn2d =  [0.75 -0.50 -0.25];
  bn2dd = [1.50 -0.50];
  bn3 =   [-0.25 0.00 0.75 0.50];
  bn3d =  [-0.75 0.00 0.75];
  bn3dd = [-1.50 0.00];
  bn4 =   [0.25 0.25 -0.25 -0.25];
  bn4d =  [0.75 0.50 -0.25];
  bn4dd = [1.50 0.50];

% Elasticity Matrix
E = (Em/((1+v)*(1-2*v)))* [(1-v)  v    v      0        0        0;
                             v  (1-v)  v      0        0        0;
                             v    v  (1-v)    0        0        0;
                             0    0    0  (1-2*v)/2    0        0;
                             0    0    0      0    (1-2*v)/2    0;
                             0    0    0      0        0   (1-2*v)/2];
  numbeamgauss=5; % Number of Gauss points for integration of beam element
  [bgpts,bgpw]=gauss(numbeamgauss);
  k_brick=zeros(24,24);% For this beam, 2 nodes, 2DOF each, is a 4 by 4
                 % matrix. 
  l=norm([x2 y2 z2]-[x1 y1 z1]);
  propertynum=num2str(element(elnum).properties);
    
  Jac=l/2;% Brick element Jacobian. 
          % Local Bending in x-y plane
  for i=1:numbeamgauss
    bricksfs=[polyval(bn1dd,bgpts(i))/Jac^2;%evaluating second
                                           %derivatives of shape
                                           %functions to use in
                                           %generating stiffness
                                           %matrix. (at gauss point)
             polyval(bn2dd,bgpts(i))/Jac;
             polyval(bn3dd,bgpts(i))/Jac^2;
             polyval(bn4dd,bgpts(i))/Jac];
             
    k_brick=k_brick+bgpw(i)*bricksfs*bricksfs'*E*Jac;%This is the Gauss
                                               %integration part. 
  end
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % 
  % Derivation of Mass matrices
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  numbeamgauss=numbeamgauss+3; %Need more gauss points for the mass
                               %matrix. 
  [bgpts,bgpw]=gauss(numbeamgauss);
  m_brick=zeros(24,24); %initialize empty mass matrix
  
  for i=1:numbeamgauss
    bricksfs=[polyval(bn1,bgpts(i));
              polyval(bn2,bgpts(i))*Jac;
              polyval(bn3,bgpts(i));
              polyval(bn4,bgpts(i))*Jac];
   
    m_brick=m_brick+bgpw(i)*bricksfs*bricksfs'*rho*Jac;% pause, and reflect
                                                 %(OK, this was for debugging)
  end
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  % Coordinate rotations
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  R1=([x2 y2 z2]-[x1 y1 z1]);% Vector along element
  lam1=R1/norm(R1);% Unit direction
  R2=([x3 y3 z3]-[x1 y1 z1]);% Unit direction to point
  R2perp=R2-dot(R2,lam1)*lam1;% Part of R2 perpendicular to lam1
  udirec=0;
  while norm(R2perp)<10*eps% If R2perp is too small, (point in line
                           % with element, we need to cover the
                           % users a$$ and generate a point that
                           % isn't. We should put out a warning,
                           % but I commented it out. 
    udirec=udirec+1;
    %disp('oops'); %This was my warning. 
    %pause
    [minval,minloc]=min(lam1);
    R2perp=zeros(1,3);
    R2perp(udirec)=1;
    R2perp=R2perp-dot(R2perp,lam1)*lam1;
  end
  %Make the unit direction vectors for rotating and put them in the
  %rotation matrix. 
  lam2=R2perp/norm(R2perp);
  lam3=cross(lam1,lam2);
  lamloc=[lam1;lam2;lam3];
  lam=sparse(12,12);
  lam(1:3,1:3)=lamloc;
  lam(4:6,4:6)=lamloc;
  lam(7:9,7:9)=lamloc;
  lam(10:12,10:12)=lamloc;

  
% $$$     lam=[lamloc z z z z z;
% $$$          z lamloc z z z z;
% $$$          z z lamloc z z z;
% $$$          z z z lamloc z z;
% $$$          z z z z lamloc z;
% $$$          z z z z z lamloc];
  element(elnum).lambda=lam;
  element(elnum).m=m;
  element(elnum).k=k;

  kg=lam'*k*lam;
  mg=lam'*m*lam;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  % Assembling matrices into global matrices
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  bn1=bnodes(1);bn2=bnodes(2);
  indices=[bn1*6+(-5:0) bn2*6+(-5:0)] ;


  K(indices,indices)=K(indices,indices)+kg;
  M(indices,indices)=M(indices,indices)+mg;

  % At this point we also know how to draw the element (what lines
  % and surfaces exist). For the beam3 element, 2 lines are
  % appropriate. Just add the pair of node numbers to the lines
  % array and that line will always be drawn.
  numlines=size(lines,1);
  lines(numlines+1,:)=[bn1 bn2];
%   lines(numlines+2,:)=[bn3 bn2];

  
  
  %If I have 4 nodes that I want to use to represent a surface, I
  %do the following.
  panelcolor=[1 0 1];% This picks a color. You can change the
                     % numbes between 0 and 1. 
  %Don't like this color? Use colorui to pick another one. Another
  %option is that if we can't see the elements separately we can
  %chunk up x*y*z, divide by x*y*x of element, see if we get
  %integer powers or not to define colors that vary by panel. 
  
  
  % You need to uncomment this line and assign values to node1,
  % node2, node3, and node4 in order to draw A SINGLE SURFACE. For
  % a brick, you need 6 lines like this. 
  %surfs=[surfs;node1 node2 node3 node4 panelcolor];
  
  %Each surface can have a different color if you like. Just change
  %the last three numbers on the row corresponding to that
  %surface. 

%diag(M)
elseif strcmp(mode,'istrainforces')
  % You don't need this
  % We need to have the stiffness matrix and the coordinate roation matrix.
 

  
elseif strcmp(mode,'draw')
elseif strcmp(mode,'buckle')
end
