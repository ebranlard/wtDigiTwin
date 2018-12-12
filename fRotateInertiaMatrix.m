function I_B = fRotateInertiaMatrix( I_A,  R_BA)
% Change coordinates system of an inertia matrix
%
% INPUTS:
%    I_A : Inertia matrix 3x3 in the coordinate system A
%    R_BA: Rotation matrix from the system A to the system B, i.e.  v|B= R_BA . v|A 
%
% AUTHOR: E.Branlard

I_B= R_BA * I_A * R_BA';
