function I_B = fTranslateInertiaMatrixFromCOG(I_G, Mass, r_BG)
% Transform inertia matrix with respect to COG to the inertia matrix with respect to point B
%
% NOTE: the vectors and the inertia matrix needs to be expressed in the same coordinate system. 
%
% INPUTS:
%    I_G  : Inertia matrix 3x3 with respect to COG
%    Mass : Mass of the body
%    r_BG: vector from point B to COG of the body
%
% AUTHOR: E.Branlard
%
if nargin==0
    I_G=eye(3); M=1;
    r_AG = [3,4,5];
    I_A  = fTranslateInertiaMatrixFromCOG (I_G, M, r_AG) ;
    I_G2 = fTranslateInertiaMatrixToCOG   (I_A, M, -r_AG);
    if norm(I_G-I_G2)>1e-8 
        fprintf('[FAIL] Inertia badly transferred to and from COG\n')
    end
    return
end

I_B= I_G  - Mass * fSkew(r_BG)^2 ;
end

function M=fSkew(x)
% returns the skew symmetric matrix M, such that: cross(x,v) = M v
M=[0 -x(3) x(2) ; x(3) 0 -x(1) ; -x(2) x(1) 0 ];
end
