function F_T_D_SBT_LH(R, a, lambda, L)
%   
% Input arguments: 
%   -a_lambda is the ratio of filament radius "a" to the wavelength
%       "lambda", i.e., a_lambda=a/lambda
%   -theta is the pitch angle
%   -nlength is the axial length of helix in the unit of its wavelength
%       for one period of helix, nlength = 1.
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
%   vphi is the mesh of helical phase, which is stored as a vector
%     
global force_LH torque_LH drag_LH
theta = atan(2*pi*R./lambda);    %Pitch angle
        
for j = 1:length(lambda)
        a_lambda(:,j) = a(:)./lambda(j);
end
        
for j = 1:length(lambda)
        nlength(j,:) = L./lambda(j);
        contour_length(j,:) = lambda(j)./cos(theta(j)) .* nlength(j,:);
end

phirange = nlength*2*pi;

res = 20;
nres=floor(nlength*res+.5);

% contour_length = lambda./cos(theta) .* nlength;
arc_length = contour_length./nres;

%Get F and T for nontranslating flagellum
omega = 1.0;    %Set rotation rate to unity
V = 0.0;        %Set velocity to zero

for i = 1:length(a)
    for j=1:length(lambda)
        for k= 1:length(L)
            [F_LH, T_LH] = SBT_helical_Lighthill(a_lambda(i,j), theta(j), nlength(j,k), res, omega, V);
            force_LH(i,j,k)=-sum(F_LH)*4*pi*arc_length(j,k);
            torque_LH(i,j,k)=sum(T_LH)*4*pi*arc_length(j,k);
        end
    end
end

%Get F and T for nonrotating flagellum
omega = 0.0;    %Set rotation rate to zero
V = 1.0;        %Set velocity to unity
for i = 1:length(a)
    for j=1:length(lambda)
        for k= 1:length(L)
            [F_LH, T_LH] = SBT_helical_Lighthill(a_lambda(i,j), theta(j), nlength(j,k), res, omega, V);
            drag_LH(i,j,k)=sum(F_LH)*4*pi*arc_length(j,k);
        end
    end
end

% save Length_Dep_SBT_Data.mat force_LH torque_LH drag_LH force_J torque_J drag_J lambda L a R