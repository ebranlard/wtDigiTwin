function [v, om, omp_v]=fComputeVelAndRotAcc(J,gxp);
    % Computes velocities at point P in the P system using the J matrix: J * GXP
    % The translation velocities are obtained from the first three lines of J
    % The rotational  velocities are obtained from the last three lines of J
    %   using the recursive formula om_dot = om1 x om2
    % J has dimensions 6 * nP
    %
    % AUTHOR: E. Branlard

    nP=size(J,2);
    if size(J,1)~=6
        error('Wrong dim for J')
    end
    v     = [0; 0; 0];
    om    = [0; 0; 0];
    omp_v = [0; 0; 0];
    for i=nP:-1:1  % NOTE Loop order is important for assessment of omp_v
        % Translation velocity of point 
        v=v+J(1:3,i)*gxp(i); 
        % gxp(J)'s contribution to OM
        om_sub=J(4:6,i)*gxp(i);
        % omp_v 
        omp_v=omp_v + cross(om_sub,om);
        % Rotational velocity of point P in P-system *)
        om=om+om_sub;
    end 
end

