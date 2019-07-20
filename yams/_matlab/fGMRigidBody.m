function MM=fGMRigidBody(Mass,J,rho)
% Generalized mass matrix for a rigid body
%
% AUTHOR: E.Branlard
    %
    S=Mass*fSkew(rho);

    MM=zeros(6,6);
    MM(1:3,1:3) = Mass*eye(3);
    MM(1:3,4:6) = -S;
    MM(4:6,1:3) = S ; % transpose(S)=-S;
    MM(4:6,4:6) = J ;
end
function M=fSkew(x)
% returns the skew symmetric matrix M, such that: cross(x,v) = M v
M=[0 -x(3) x(2) ; x(3) 0 -x(1) ; -x(2) x(1) 0 ];
end
