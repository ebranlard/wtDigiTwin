function I_B = fTranslateInertiaMatrix(I_A, Mass, r_BG, r_AG)
% Transform inertia matrix with respect to point A to the inertia matrix with respect to point B
%
% NOTE: the vectors and the inertia matrix needs to be expressed in the same coordinate system. 
% NOTE: one of the vector r_BG or r_AG may be empty or 0 instead of [0,0,0];
% NOTE: if r_AG is not provided it is assumed to be 0, i.e. A=G
% To avoid this confusion you can use fTranslateInertiaMatrixFromCOG  and fTranslateInertiaMatrixToCOG
%
% INPUTS:
%    I_A  : Inertia matrix 3x3 in the coordinate system A
%    Mass : Mass of the body
%    r_BG: vector from point B to COG of the body
% 
% OPTIONAL INPUTS:
%    r_AG: vector from point A to point G
%
% AUTHOR: E.Branlard
%
if nargin==0
    %% Simple transfer along perpendicular axis
    I_G=eye(3);
    M=1;
    r_OG = [ 1; 2; 10 ];
    r_OA = r_OG + [0; 0 ; 2];
    r_OB = r_OG + [0; 0 ;-3];
    r_AG = r_OG-r_OA;
    r_BG = r_OG-r_OB;
    I_A  = fTranslateInertiaMatrix (I_G, M, r_AG)    ;
    I_B  = fTranslateInertiaMatrix (I_G, M, r_BG)    ;
    I_G2 = fTranslateInertiaMatrix (I_A, M, 0, -r_AG);
    I_G3 = fTranslateInertiaMatrix (I_B, M, 0, -r_BG);
    I_B2 = fTranslateInertiaMatrix (I_A, M, r_BG, r_AG);
    I_A2 = fTranslateInertiaMatrix (I_B, M, r_AG, r_BG);
    % Performing tests 
    if abs(I_A(1,1)- ( I_G(1,1)  + M * 2^2))>1e-8  
        fprintf('[FAIL] Inertia badly transferred from COG to point A\n')
    end
    if abs(I_B(1,1)- ( I_G(1,1)  + M * 3^2))>1e-8  
        fprintf('[FAIL] Inertia badly transferred from COG to point B\n')
    end
    if norm(I_A-I_A2)>1e-8 && norm(I_B-I_B2)>1e-8
        fprintf('[FAIL] Inertia badly transferred from point A to point B\n')
    end
    if norm(I_G-I_G2)>1e-8 && norm(I_G-I_G3)>1e-8
        fprintf('[FAIL] Inertia badly transferred from points to COG\n')
    end
    %% Generic transfer
    I_G=diag([1,2,3]); M=1;
    r_OG = [ 1; 2; 10 ];
    r_OA = r_OG + [5;  8 ; 2];
    r_OB = r_OG + [4; -6 ;-3];
    r_AG = r_OG-r_OA;
    r_BG = r_OG-r_OB;
    I_A  = fTranslateInertiaMatrix (I_G, M, r_AG)      ;
    I_B  = fTranslateInertiaMatrix (I_G, M, r_BG)      ;
    I_B2 = fTranslateInertiaMatrix (I_A, M, r_BG, r_AG);
    I_A2 = fTranslateInertiaMatrix (I_B, M, r_AG, r_BG);
    % Performing tests 
    if norm(I_A-I_A2)>1e-8 && norm(I_B-I_B2)>1e-8
        fprintf('[FAIL] Inertia badly transferred from point A to point B\n')
    end
    return
end



% Setting up default values
if ~exist('r_AG','var'); r_AG=[0;0;0]; end; % We assume A=G !!!
% Allowing the user to give empty or 0 inputs for a zero vector
if length(r_AG)<3; r_AG=[0;0;0]; end;
if length(r_BG)<3; r_BG=[0;0;0]; end;

I_B= I_A  - Mass * ( fSkew(r_BG)^2  - fSkew(r_AG)^2 )  ;

end

function M=fSkew(x)
% returns the skew symmetric matrix M, such that: cross(x,v) = M v
M=[0 -x(3) x(2) ; x(3) 0 -x(1) ; -x(2) x(1) 0 ];
end
