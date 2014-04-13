function [F, T] = SBT_helical_Lighthill(a_lambda, theta, nlength, res, omega, V) 
%   
% Input arguments: 
%   -a_lambda is the ratio of filament radius "a" to the wavelength
%       "lambda", i.e., a_lambda=a/lambda
%   -theta is the pitch angle
%   -nlength is the axial length of helix in the units of wavelength
%       for one period of the helix, nlength = 1.
%       when nlength is a vector, it stores the length of a set of helixes
%       to compute for the force. 
%   -res is the number of meshes per wave length, typally I set res=20.
%    
%   -omega is the linear rotation speed
%   -V is the translation speed
%   
% Output arguments:
%   -F, T are the mean forces per arclength along the translation and rotation
%       directions, respectively. which are stored as vectors of the same size as
%       as nlength.
%       *Note that these output forces are normalized by 4 pi mu R
%       
    

    function m = rotmatrix(angle)
        m=[cos(angle), sin(angle), 0; -sin(angle), cos(angle), 0; 0, 0, 1];
    end
F=zeros(size(nlength));
T=zeros(size(nlength)); 

for k=1:numel(nlength)
    
phimin=0; phimax=nlength(k)*(2*pi);
xi=@(phi) sqrt(2 - 2*cos(phi(2)-phi(1))+(phi(2)-phi(1))^2/(tan(theta)^2));
fx12=@(phi) [cos(phi(1))-cos(phi(2)), sin(phi(1))-sin(phi(2)), (phi(1)-phi(2))/tan(theta)];
px12=@(phi) (rotmatrix(phi(1))*(fx12(phi)'*fx12(phi))*rotmatrix(-phi(2))/xi(phi)^3+rotmatrix(phi(1)-phi(2))/xi(phi));
nres=floor(nlength(k)*res+.5);
dphi=(phimax-phimin)/nres;
vphi=phimin+.5*dphi:dphi:phimax-.5*dphi;
mattr=zeros(3*nres, 3*nres);
cutofflen=a_lambda*pi*exp(0.5)*cos(theta);
cl=log(.5*dphi/cutofflen);
for i=1:nres
    
mattr(i,i)=1+cl;
mattr([nres+i, 2*nres+i], [nres+i, 2*nres+i])=[cos(theta)^2+(1+sin(theta)^2)*cl, (cl-1)*sin(theta)*cos(theta); (cl-1)*cos(theta)*sin(theta), sin(theta)^2+(1+cos(theta)^2)*cl];
end

for i=1:nres
for j=1:nres
if (abs(vphi(i)-vphi(j))<cutofflen+.5*dphi)
continue;
end
mattr([i, nres+i, 2*nres+i], [j, nres+j, 2*nres+j])=px12([vphi(i), vphi(j)])*dphi*.5/sin(theta);
end
end
v=zeros(1, 3*nres);
v(nres+1:2*nres)=omega;
v(2*nres+1:3*nres)=V;
f=mattr\v';
F(k)=sum(f(2*nres+1:3*nres));%/nres;
T(k)=sum(f(nres+1: 2*nres));%/nres;
end
end

