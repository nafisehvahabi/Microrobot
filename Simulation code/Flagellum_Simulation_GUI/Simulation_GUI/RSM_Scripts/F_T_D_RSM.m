function F_T_D_RSM(R, a, lambda, L)

% Input Arguments
%   -R the helical radius (single valued)
%   -a the filament radius (vector)
%   -lambda helical pitch or wavelength (vector)
%   -L axial length of the helix (vector)

global force_RSM torque_RSM drag_RSM

zz = 250;

for i = 1:length(a)
    for j=1:length(lambda)
        for k= 1:length(L)
            
           [F1,F2,F3,T1,T2,T3,A1,A2,A3] = ExpFlagellumFullTubeMatrix(zz,a(i),R,L(k),lambda(j));
            force_RSM(i,j,k) = F1;
            force_RSM_Y(i,j,k) = F2;  %% Y and Z values are not used
            force_RSM_Z(i,j,k) = F3;

            torque_RSM(i,j,k) = T1;
            torque_RSM_Y(i,j,k) = T2;
            torque_RSM_Z(i,j,k) = T3;

            drag_RSM(i,j,k) = A1;
            drag_RSM_Y(i,j,k) = A2;
            drag_RSM_Z(i,j,k) = A3;
          
        end
    end
end


