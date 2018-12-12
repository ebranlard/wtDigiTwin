classdef Connection < handle
%% Documentation
% Class to handle the connection between Bodies
% 
% AUTHOR: E. Branlard

properties
    Name;
    Type;
    s_C_0_inB; % Initial position of connection point in local coord
    s_C_inB  ; % Position of connection point in local coord
    s_C      ; % Position of connection point in global coord
    R_ci_0   ; % Initial orientation
    R_ci     ; % Rotation from connection point to chilld body
    B_ci=[];
    nj; % number of DOFs required
    I_DOF; % Index of DOFs in global vector
    ParentNode;
    % --- Spherical
    JointRotations;
%     r_C     ; % Position of connection point in global coordinates
% ac.connectTo(Sft,'Type','SphericalJoint','Point_inB',r_NS_inN,'JointRotations',{'z'});
end
methods
    function o=Connection(varargin)
        p = fInputParser();
        % --- Key, value parameters
        p.addParameter('Name'          ,'',@ischar   );
        p.addParameter('Type'          ,'',@ischar   );
        p.addParameter('JointRotations','',@iscell   );
        p.addParameter('RelPoint'      ,[],@isnumeric);
        p.addParameter('ParentNode'    ,[],@isnumeric);
        p.addParameter('Orientation'   ,eye(3),@isnumeric);
        p.parse(varargin{:});
        p=p.Results;
        o.Name      = p.Name       ;
        o.Type      = p.Type       ;
        o.s_C_0_inB = p.RelPoint   ;
        o.s_C_inB   = o.s_C_0_inB    ;
        o.R_ci_0    = p.Orientation;
        o.R_ci      = o.R_ci_0     ;
        o.ParentNode=p.ParentNode;
        switch o.Type
            case 'Rigid'
                o.nj=0;

            case 'SphericalJoint'
                o.JointRotations=p.JointRotations;
                o.nj=length(o.JointRotations);
            otherwise 
                error('Unknown connection type %s',o.Type)
        end
    end

    function updateKinematics(j,q,qdot)
        myq    = q   (j.I_DOF);
        myqdot = qdot(j.I_DOF);

        j.B_ci=zeros(6,j.nj);

        switch j.Type
            case 'Rigid'
                %j.R_ci=j.R_ci_0;
            case 'SphericalJoint'
                R=eye(3);
                for ir=1:length(j.JointRotations)
                    switch j.JointRotations{ir}
                        case 'x'; I=[1;0;0];
                        case 'y'; I=[0;1;0];
                        case 'z'; I=[0;0;1];
                        otherwise; error('Here');
                    end
                    % Setting Bhat column by column
                    j.B_ci(4:6,ir) = R*I; % NOTE: needs to be done before R updates
                    % Updating rotation matrix
                    R=R*fRot(j.JointRotations{ir}, myq(ir));
                end
                j.R_ci=R*j.R_ci_0;

                
            otherwise 
                error('Unknown connection type %s',o.Type)
        end
    end



end % methods

end % class
