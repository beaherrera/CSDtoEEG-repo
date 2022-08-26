function [Ax,Ay]=calRotMatriz(Nor)
% author: Dr. Jorge Riera

eps=1e-9;
if Nor(1) >= 0
 psi=atan2(-Nor(2),Nor(3));  
 Ax=[1,0,0;0,cos(psi),sin(psi);0,-sin(psi),cos(psi)];
 N_x=Ax*Nor';
 cita=acos(N_x(3)/sqrt(Nor(1)^2+N_x(3)^2));
 Ay=[cos(cita),0,-sin(cita);0,1,0;sin(cita),0,cos(cita)];
 N_xy=Ay*N_x;
else
 psi=atan2(Nor(2),Nor(3));   
 Ax=[1,0,0;0,cos(psi),-sin(psi);0,sin(psi),cos(psi)];
 N_x=Ax*Nor';
 cita=acos(N_x(3)/sqrt(Nor(1)^2+N_x(3)^2));
 Ay=[cos(cita),0,sin(cita);0,1,0;-sin(cita),0,cos(cita)];
 N_xy=Ay*N_x;  
end

if ((abs(N_xy(1)) > eps) || (abs(N_xy(2)) > eps)), disp('Error Calculating Rotation Matrices'); end

end