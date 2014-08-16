function BuildMatrixRegStokes3D(x,X0,N,N0,d)
%%
%% BUILD THE MATRIX A FOR STOKES FLOW
%%

% d is the regularization parameter
% miu is the viscosity

%  the Kernal Matrix as the global variable
global A

mu = 1;
fac0 = 1/(8*pi*mu);
%fac0 = 1;

AXX = zeros(N,1);
AYY = zeros(N,1);
AZZ = zeros(N,1);

AXY = zeros(N,1);
AXZ = zeros(N,1);
AYZ = zeros(N,1);

N3 = 3*N;
d2 = d^2;

for i=1:N0
  dx = X0(i,1)-x(:,1); dy = X0(i,2)-x(:,2); 
  dz = X0(i,3)-x(:,3);
  %dz = dx - dx;
  r2 = dx.^2 + dy.^2 + dz.^2;

%% FOR STOKES FLOWS
 tmp = (r2+d2).^(3/2);
 H2 = 1./tmp;
 H1 = (r2+2*d2).*H2;
 

%--------------------------

  AXX = H1 + H2.*dx.*dx ;
  AYY = H1 + H2.*dy.*dy ;
  AZZ = H1 + H2.*dz.*dz;
  AXY = H2.*dx.*dy;
  AXZ = H2.*dx.*dz;
  AYZ = H2.*dy.*dz;

  A(3*i-2,1:3:N3) = AXX;
  A(3*i-2,2:3:N3) = AXY;
  A(3*i-2,3:3:N3) = AXZ;
  
  A(3*i-1,1:3:N3) = AXY;
  A(3*i-1,2:3:N3) = AYY;
  A(3*i-1,3:3:N3) = AYZ;
  
  A(3*i,1:3:N3) = AXZ;
  A(3*i,2:3:N3) = AYZ;
  A(3*i,3:3:N3) = AZZ;
  
end
A = A*fac0;
end



