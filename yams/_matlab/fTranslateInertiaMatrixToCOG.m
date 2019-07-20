function I_G = fTranslateInertiaMatrixToCOG(I_B, Mass, r_GB)
% Transform inertia matrix with respect to point B to the inertia matrix with respect to the COG
%
% NOTE: the vectors and the inertia matrix needs to be expressed in the same coordinate system. 
%
% INPUTS:
%    I_G  : Inertia matrix 3x3 with respect to COG
%    Mass : Mass of the body
%    r_GB: vector from COG to point B 
%
% AUTHOR: E.Branlard
if nargin==0
    I_A=eye(3); M=1;
    r_GA = [3,4,5];
    I_G = fTranslateInertiaMatrixToCOG   (I_A, M, r_GA);
    I_A2= fTranslateInertiaMatrixFromCOG (I_G, M,-r_GA) ;
    if norm(I_A-I_A2)>1e-8 
        fprintf('[FAIL] Inertia badly transferred to and from COG\n')
    end
    return
end


I_G= I_B  + Mass * fSkew(r_GB)^2 ;
end

function M=fSkew(x)
% returns the skew symmetric matrix M, such that: cross(x,v) = M v
M=[0 -x(3) x(2) ; x(3) 0 -x(1) ; -x(2) x(1) 0 ];
end
