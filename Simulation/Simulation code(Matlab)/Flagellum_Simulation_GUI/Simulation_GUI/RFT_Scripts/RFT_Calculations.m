function[force_GH_RFT, torque_GH_RFT, drag_GH_RFT,force_LH_RFT torque_LH_RFT drag_LH_RFT] = RFT_Calculations(R, a, lambda, L)

%Gray and Hancock resistive force theory
%   
% Input arguments: 
%   -R helical radius
%   -a filament radius
%   -lambda is the helical pitch
%   -L is axial length
%   
% Output arguments:
%   -force_GH_RFT torque_GH_RFT drag_GH_RFT are the Gray and Hancock RFT axial force, torque and drag, respectively
%   -force_LH_RFT torque_LH_RFT drag_LH_RFT are the Lighthill RFT axial force, torque and drag, respectively
%       *Note that these output forces are normalized by 4 pi mu R
%     
% global force_GH_RFT torque_GH_RFT drag_GH_RFT
% global force_LH_RFT torque_LH_RFT drag_LH_RFT

theta=atan(2*pi./lambda);  % Determine pitch angle

for i = 1:length(a)
    for j=1:length(lambda)
        for k= 1:length(L)
                
            %%  Gray and Hancock RFT
            Ct=2*pi/[log(2*lambda(j)/a(i))-.5];  %No mu because will be nondimensionalized
            Cn=4*pi/[log(2*lambda(j)/a(i))+.5];

            force_GH_RFT(i,j,k)=(Cn-Ct)*sin(theta(j))*L(k); %No omega*R
            torque_GH_RFT(i,j,k)=[Cn*cos(theta(j))^2+Ct*sin(theta(j))^2]*L(k)/cos(theta(j)); %No omega*R^2
            drag_GH_RFT(i,j,k)=[Cn*sin(theta(j))^2+Ct*cos(theta(j))^2]*L(k)/cos(theta(j)); %No omega*R^2

            %%  Lighthill RFT
            Ct=2*pi/[log(0.18*lambda(j)/a(i)/cos(theta(j)))];  %No mu because nondimensionalized
            Cn=4*pi/[log(0.18*lambda(j)/a(i)/cos(theta(j)))+.5];

            force_LH_RFT(i,j,k)=(Cn-Ct)*sin(theta(j))*L(k); %No omega*R
            torque_LH_RFT(i,j,k)=[Cn*cos(theta(j))^2+Ct*sin(theta(j))^2]*L(k)/cos(theta(j)); %No omega*R^2
            drag_LH_RFT(i,j,k)=[Cn*sin(theta(j))^2+Ct*cos(theta(j))^2]*L(k)/cos(theta(j)); %No U*R
        end
    end
end

