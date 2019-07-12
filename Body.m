%
classdef Body < handle
% classdef Body < matlab.mixin.Copyable;

%% Documentation

properties
    Name;
    Type;
    %
    R_0b;  % Transformation matrix from body to ground
    R_pb;  % Transformation matrix from body to ground
    RR_0b; % Transformation matrix from body to ground for vectors of size 6 + nf
    r_O; % Position of Origin in global coordinates
    % Local coordinate "primed variables"
    v_O_inB;      % Linear velocity of origin in local coordinate system 
    om_O_inB;     % Angular velocity of origin in local coordinate system 
    a_O_inB;      % Linear acceleration of origin in local coordinate system 
    a_O_v_inB;    % Linear acceleration due to "velocity DOF"
    a_O_a_inB;    % Linear acceleration due to "acceleration DOF"
    omp_O_inB;    % Angular acceleration of origin in local coordinate system 
    omp_O_v_inB;  % Angular acceleration due to "velocity DOF"
    omp_O_a_inB;  % Angular acceleration due to "acceleration DOF"
    s_G_inB;  %

    % --- Connection point
    Children;
    Connections;
    %
    B;
    B_inB;
    BB_inB;
    I_DOF; % Index of DOFs in global vector
%     r_C     ; % Position of connection point in global coordinates
%     s_C_inB ; % Position of connection point in local coord
%     s_C0_inB; % Initial position of connection point in local coord

    %  Environment
    Mtop;
    gravity;


    % --- Rigid body
    Mass;
    J_G_inB;
    J_O_inB;

    % --- Beam Flexible body related
    s_span; % SpanLine for shape functions
    nSpan;
    % Material/bending/inertia properties
    GKt; 
    EIy; % Flexural rigidity /bending stiffness [Nm^2]
    EIz; % Flexural rigidity /bending stiffness [Nm^2]
    m    ; % Mass Distribution as function of s_span
    jxxG; %  [Kg m]
    s_P;
    s_P0;
    s_G;
    s_G0;
    rho_G;      % Position of COG in Beam coordinate
    rho_G0;     % Initial Position of COG in Beam coordinate
    rho_G0_inS; % Initial Position of COG in cross section
    %
    nf; % Number of shape functions
    PhiU; % Shape function deflection @mean line
    PhiV; % Shape function slope      @mean line
    PhiK; % Shape function curvature  @mean line
    PhiUG;% Shape function deflection @cog
    V0; % Undeflected slope   %NOTE U0=s_P0
    K0; % Undeflected curvature  
    ModeOmegas;
    gzf   ; % Generalized Coordinates for flexible motion
    gzpf  ; % Generalized Coordinates for flexible motion
    gzppf ; % Generalized Coordinates for flexible motion
    U; % Deflected shape of mean axis
    V;
    K;
    V_tot; % V+V0
    K_tot; % K+K0
    UP; % Deflected velocity of mean axis
    P_inertia; % Loads along the beam due to inertia
    P_ext;     % Loads along the beam due to external forces
    P_tot; % total loading
    Fx;
    bOrthogonal;
    bAxialCorr;
    bStiffnessFromGM;
    bStiffeningMtop;
    bStiffeningSelfWeight;
    DampingMethod='None';
    DampingParams;

    %
    % --- Shape Integrals
    sigma     ; % I1 ; % 3x1   -
    Psi       ; % S  ; % 3xnf  -
    Upsilon_kl; % Skl; % nfxnf -
    Sigma_kl  ; % Ikl; % nfxnf -
    sigma_kl  ; % ikl; % nfxnf -
    I_Jxx;
    GM_Jxx;
    ImomX;
    M1;
    
    % Generalized coordinates variables -  Matrices
    MM;
    KK;
    KKg;  % Geometrical stiffening
    KK0;
    DD;
    GF;

    % Shape Integrals 
    bUseShapeIntegrals;
    Smom;
    Imom;
    Smom1;
    IW;    % Span Integration Weight  
    IW_x;  % Span Integration Weight
    IW_m;  % Mass Integration Weight
    IW_xm; % Mass Integration Weight
    IW_U;  % Shape function Integration Weight 

end
properties(SetAccess = private, Hidden = true)
    bCompatibility;
end
% properties(SetAccess = protected, Hidden = true)
% end
% properties(Abstract = true, Hidden=true);
% end

% --------------------------------------------------------------------------------
% --- Main public functions: constructor, read and write
% --------------------------------------------------------------------------------
methods
    function o=Body(varargin)
        %%  Setting default values
        o.setDefault();
        %% Optional arguments as Key-values
        p = fInputParser();
        % --- Key, value parameters
        p.addParameter('Name','',@ischar);
        p.addParameter('Type','',@ischar);
        p.parse(varargin{:});
        p=p.Results;
        %% Hack to run unit tests
        if nargin==0
            [~,ParentDir]=fileparts(pwd);
            if isequal('YAMS',ParentDir)
                run Tests;
            end
        else
            o.Name=p.Name;
            o.Type=p.Type;
        end

    end
    function o=setDefault(o)
        o.bCompatibility=false;
        o.bUseShapeIntegrals=true;
        o.bAxialCorr=true;
        o.bStiffnessFromGM=false;
        o.bStiffeningMtop       = false;
        o.bStiffeningSelfWeight = false;
        o.ModeOmegas=[];
        o.Mtop=0;
        o.gravity=9.81;
    end

    function o=init(o)
        switch o.Type
            case 'Ground'
                o.nf=0;
                updateKinematics(o,[0;0;0],eye(3),zeros(o.nf,1),zeros(6+o.nf,1),zeros(6+o.nf,1)); % IMPORTANT HERE, R_0b=eye(3)
            case 'Rigid'
                % Checks
                o.check_not_empty('Mass','s_G_inB');
                if isempty(o.J_G_inB) && isempty(o.J_O_inB); error('No inertia provided'); end;
                if isempty(o.J_G_inB)
                    o.J_G_inB=fTranslateInertiaMatrixToCOG(o.J_O_inB, o.Mass, -o.s_G_inB);
                end
                if isempty(o.J_O_inB)
                    o.J_O_inB=fTranslateInertiaMatrixFromCOG(o.J_G_inB, o.Mass, o.s_G_inB);
                end
                o.nf=0;
                gz0=zeros(o.nf,1);
                o.bUseShapeIntegrals=false;
                o.bAxialCorr=false;
            case 'FlexibleBeam'
                % INPUTS REQUIRED FOR FLEXIBLE BODY
                %o.check_not_empty('s_span','m','PhiU');
                o.check_not_empty('s_span','m')
                if isempty(o.EIy) && isempty(o.EIz); error('Provide at least EIy or EIZ'); end;
                % Useful variables
                o.nf = length(o.PhiU);
                o.nSpan=length(o.s_span);
                % --- Optional arguments
                if isempty(o.PhiUG); o.PhiUG = o.PhiU; end;
                if isempty(o.gzf); o.gzf=zeros(o.nf,1); end;
                if isempty(o.EIz ); o.EIz  = zeros(1,o.nSpan); end;
                if isempty(o.EIy ); o.EIy  = zeros(1,o.nSpan); end;
                if isempty(o.GKt ); o.GKt  = zeros(1,o.nSpan); end;
                if isempty(o.jxxG); o.jxxG = zeros(1,o.nSpan); end;
                if isempty(o.s_P0) % Points P0 - Undeformed points of the body
                    o.s_P0=zeros(3,o.nSpan);
                    o.s_P0(1,:)=o.s_span; 
                    o.s_P0(2,:)=0;
                    o.s_P0(3,:)=0;
                end
                % If V&K not provided, num computation, otherwise check their values
                if o.nf>0
                    [o.PhiV,o.PhiK] = fBeamSlopeCurvature(o.s_span,o.PhiU,o.PhiV,o.PhiK,1e-2);
                    [o.V0,o.K0]     = fBeamSlopeCurvature(o.s_span,o.s_P0,o.V0,o.K0,1e-2)    ;
                    if isempty(o.s_G0); o.s_G0=o.s_P0; end;
                    if isempty(o.rho_G0_inS); o.rho_G0_inS=zeros(3,o.nSpan); end;
                    if isempty(o.rho_G0    ); 
                        o.rho_G0 =zeros(3,o.nSpan);
                        for i=1:o.nSpan
                            o.rho_G0(1:3,i) =fRotx(o.V0(1,i))*o.rho_G0_inS(:,i);
                        end;
                    end
                end
                % --- Sanity
                o.line_vector_me('s_span','m','EIy','EIz','GKt','jxxG');
                %o.check_same_size('s_span','m','EIz');

                % --- Initial conditions 
                gz0=o.gzf;
            otherwise
                error('Unknown body type %s',o.Type);
        end

        if ~isequal(o.Type,'Ground')
            % Setting initial position (x_0,R_0b,gz,v_0,a_v_0)
             %  updateKinematics(o,x_0,R_0b,gz,v_0,a_v_0)
            updateKinematics(o,[0;0;0],zeros(3,3),gz0,zeros(6+o.nf,1),zeros(6+o.nf,1));

            % Integration weights are always useful
            o.computeIntegrationWeights();

            if o.bUseShapeIntegrals
                o.computeShapeIntegrals();
            end

            % Matrices
            o.computeMassMatrix();
            o.computeStiffnessMatrix();
            o.computeDampingMatrix();
        end

    end

    function updateKinematics(o,x_0,R_0b,gz,v_0,a_v_0)
        % Updating position of body origin in global coordinates
        o.r_O = x_0(1:3);
        o.gzf = gz(:);
        % Updating Transformation matrix
        o.R_0b=R_0b;
        o.RR_0b = blkdiag(o.R_0b, o.R_0b , eye(o.nf));
        % Updating velocity
        v_inB         = (o.RR_0b')*v_0  ;
        o.v_O_inB     = v_inB(1:3);
        o.om_O_inB    = v_inB(4:6);
        o.gzpf        = v_inB(6+(1:o.nf));
        % Updating acceleration
        a_v_inB       = (o.RR_0b')*a_v_0; 
        o.a_O_v_inB   = a_v_inB(1:3)  ;
        o.omp_O_v_inB = a_v_inB(4:6)  ;
        o.gzppf       = a_v_inB(6+(1:o.nf))  ;

        % --- Calculation of deformations wrt straight beam axis, curvature (K) and velocities (UP)
        if isequal(o.Type,'FlexibleBeam') && o.nf>0
            % Deflections shape
            o.U  = zeros(3,o.nSpan);
            o.V  = zeros(3,o.nSpan);
            o.K  = zeros(3,o.nSpan);
            %o.U(1,:) = o.s_span; 
            o.UP = zeros(3,o.nSpan);
            for j=1:o.nf
                o.U (1:3,:) = o.U(1:3,:)  + o.gzf(j)  * o.PhiU{j}(1:3,:);
                o.UP(1:3,:) = o.UP(1:3,:) + o.gzpf(j) * o.PhiU{j}(1:3,:);
                o.V (1:3,:) = o.V(1:3,:)  + o.gzf(j)  * o.PhiV{j}(1:3,:);
                o.K (1:3,:) = o.K(1:3,:)  + o.gzf(j)  * o.PhiK{j}(1:3,:);
            end
            o.V_tot=o.V+o.V0;
            o.K_tot=o.K+o.K0;

            % Poisiton of mean line
            o.s_P=o.s_P0+o.U;

            % Position of deflected COG
            o.rho_G      = zeros(3,o.nSpan);
            o.rho_G(2,:) = o.rho_G0_inS(2,:).*cos(o.V_tot(1,:))-o.rho_G0_inS(3,:).*sin(o.V_tot(1,:));
            o.rho_G(3,:) = o.rho_G0_inS(2,:).*sin(o.V_tot(1,:))+o.rho_G0_inS(3,:).*cos(o.V_tot(1,:));
            o.s_G = o.s_P+o.rho_G;
            % Alternative:
            %rho_G2     = zeros(3,o.nSpan);
            %rho_G2(2,:) = o.rho_G0(2,:).*cos(o.V(1,:))-o.rho_G0(3,:).*sin(o.V(1,:));
            %rho_G2(3,:) = o.rho_G0(2,:).*sin(o.V(1,:))+o.rho_G0(3,:).*cos(o.V(1,:));
            %compare(o.rho_G,rho_G2,'rho_G');

            % Position of connection point
            for ic=1:length(o.Connections)
                iNode=o.Connections{ic}.ParentNode;
                %o.Connections{ic}.s_C_inB = o.U(1:3,iNode);
                o.Connections{ic}.s_C_inB = o.s_P(1:3,iNode);
            end
        else
            % No update of connection point for rigid bodies
        end
     end % updateKinematics

    function o=computeInertiaForces(o)
        om  = o.om_O_inB ;
        a   = o.a_O_v_inB; % NOTE : THAT'S A GTILDE !!!!
        omp = o.omp_O_v_inB;

        % --- Computing inertia term - All vectors are expressed in the body frame!
        % Acceleration : om x (om x s) + 2 om x v + omp x s
        A=zeros(3,o.nSpan);
        s1 = o.s_P(1,:);  % Position of point in body frame
        s2 = o.s_P(2,:);  % Position of point in body frame
        s3 = o.s_P(3,:);  % Position of point in body frame
        v1 = o.UP(1,:); % Velocity of point in body frame 
        v2 = o.UP(2,:); % Velocity of point in body frame 
        v3 = o.UP(3,:); % Velocity of point in body frame 
        A(1,:) = a(1) - ( om(2)*(om(1)*s2-s1*om(2))-om(3)*(om(3)*s1-om(1)*s3) + 2*(om(2)*v3-v2*om(3))  +omp(2)*s3-omp(3)*s2);
        A(2,:) = a(2) - ( om(3)*(om(2)*s3-s2*om(3))-om(1)*(om(1)*s2-om(2)*s1) + 2*(om(3)*v1-v3*om(1))  +omp(3)*s1-omp(1)*s3);
        A(3,:) = a(3) - ( om(1)*(om(3)*s1-s3*om(1))-om(2)*(om(2)*s3-om(3)*s2) + 2*(om(1)*v2-v1*om(2))  +omp(1)*s2-omp(2)*s1);

        % Total inertia loads
        o.P_inertia=zeros(3,o.nSpan);
        o.P_inertia(1,:) = o.m .* A(1,:);
        o.P_inertia(2,:) = o.m .* A(2,:);
        o.P_inertia(3,:) = o.m .* A(3,:);

        % Integration of axial loading Fx(r) = int_r^R m(r).ax(r) dr
        o.Fx= fcumtrapzlr(o.s_span, o.P_inertia(1,:)) ;

        % --- Axial load from acceleration
        % Method 1
        Corr=zeros(3,o.nSpan);
        Corr(2,:) = o.Fx.*o.K(2,:)  - o.P_inertia(1,:).*o.V(2,:);
        Corr(3,:) = o.Fx.*o.K(3,:)  - o.P_inertia(1,:).*o.V(3,:);

        % Method 2
        KKCorr=zeros(o.nf,o.nf);
        for j=1:o.nf
%             xx=trapz(o.s_span, o.PhiV{i}(1,:).* o.PhiV{j}(1,:));
            yy=trapz(o.s_span, o.PhiV{j}(2,:).* o.PhiV{j}(2,:).*o.Fx);
%             zz=trapz(o.s_span, o.PhiV{i}(3,:).* o.PhiV{j}(3,:));
            KKCorr(j,j)=yy;
        end
        if any(isnan(KKCorr))
            disp('Probelm in KKcorr computation')
            keyboard
        end
        o.KK(7:end,7:end)= o.KK0(7:end,7:end)+KKCorr*0;

        o.P_inertia = o.P_inertia + Corr;
    end


    function computeGeneralizedForces(o)
      % Calculation of generalized forces in each DOF 

      % Integration of axial loading
      % Fx(r) = int_r^R m(r).ax(r) dr
      o.Fx= fcumtrapzlr(o.s_span, o.P_inertia(1,:)) ;

      % Total force % TODO add external forces
      o.P_tot = o.P_inertia;

      o.GF=zeros(1,6+o.nf);
      o.GF(1)=o.Fx(1);
      for I=1:o.nSpan
          PK=o.P_tot(2,I);
          PF=o.P_tot(3,I);
          % Forces
          o.GF(2)=o.GF(2)+o.IW(I)*PK;
          o.GF(3)=o.GF(3)+o.IW(I)*PF;
          % Moments
          %GF(4)=// Add torsional moment
          o.GF(5)=o.GF(5)-o.IW_x(1,I)*PF;
          o.GF(6)=o.GF(6)+o.IW_x(1,I)*PK;
      end;
      %(* korrektion for usymmetri af tvrkraft fra vinkelacceleration *)
      %GFV[2]:=GFV[2]-(Smom-Smom1)*OMPZ;   (* Fy *)                       (* 18.9.94 *)
      %GFV[3]:=GFV[3]+(Smom-Smom1)*OMPY;   (* Fz *)
      % Elastic DOFs -  \int \u . \p ds 
      if ~isempty(o.IW_U)
          for j=1:o.nf; o.GF(6+j)=sum(            dot(o.IW_U{j} , o.P_tot)); end; 
      else 
          for j=1:o.nf; o.GF(6+j)=trapz(o.s_span, dot(o.PhiU{j}, o.P_tot)); end;
      end
    end


    function o=computeIntegrationWeights(o)
        if isequal(o.Type,'FlexibleBeam')
            % Integration weight factors
            [o.IW,o.IW_x,o.IW_m,o.IW_xm,o.IW_U,o.Mass,o.Smom,o.Smom1] = fIntegrationWeightsBeam(o.s_span,o.m,o.PhiU,o.PhiV,o.PhiK);
        end
    end

    function o=computeShapeIntegrals(o)
        % Compute the shape integral
        if ~isempty(o.sigma); warning('Recomputing shape integrals. Not necessary'); end
        if o.nf>0
            if isempty(o.PhiV); error('PhiV shoud have been set by procedure init. Did you call it?'); end
            if isempty(o.PhiK); error('PhiK shoud have been set by procedure init. Did you call it?'); end
        end

        % NOTE: using IW_xm will change the results a bit for Mtt and Mtg
        [o.Psi,o.Upsilon_kl,o.Sigma_kl,o.sigma_kl,o.sigma,o.Mass,o.ImomX,o.I_Jxx,o.GM_Jxx,o.M1] = ...
            fShapeIntegralsBeam(o.s_G0,o.s_span,o.m,o.jxxG,o.PhiUG,o.PhiV,'IW_xm',o.IW_xm,'bAxialCorr',o.bAxialCorr,'V0',o.V0);
            %fShapeIntegralsBeam(o.s_G0,o.s_span,o.m,o.jxxG,o.PhiUG,o.PhiV,'bAxialCorr',o.bAxialCorr,'V0',o.V0);
    end

    function o=computeMassMatrix(o)
        switch o.Type
            case 'Rigid'
                % NOTE: function of gzf!
                o.MM = fGMRigidBody(o.Mass,o.J_O_inB,o.s_G_inB);
            case 'FlexibleBeam'
                if o.bUseShapeIntegrals
                    o.MM = fGMBeamShapeIntegrals(o.gzf,o.Psi,o.Upsilon_kl,o.Sigma_kl,o.sigma_kl,o.sigma,o.Mass,o.ImomX,o.I_Jxx,o.GM_Jxx,o.bOrthogonal,o.bAxialCorr,o.M1);
                else
                    o.MM = fGMBeamStandalone     (o.s_G  ,o.s_span,o.m,o.jxxG,o.PhiUG,o.PhiV,o.IW,o.IW_xm,o.bOrthogonal,o.bAxialCorr,o.V_tot);
                end
        end
    end

    function o=computeStiffnessMatrix(o)
        switch o.Type
            case 'Rigid'
                o.KK=zeros(6,6);
            case 'FlexibleBeam'
            %[o.KK0,o.KKg] = fGKBeam(o.s_span,o.EI*0,o.EI,o.EI,o.PhiK,o.bOrthogonal);
            [o.KK0,o.KKg] = fGKBeam(o.s_span,o.GKt,o.EIy,o.EIz,o.PhiK,o.PhiV,o.V0,o.K0,'bOrth',o.bOrthogonal,'bStiffeningMtop',o.bStiffeningMtop,'bStiffeningSelfWeight',o.bStiffeningSelfWeight,'Mtop',o.Mtop,'m',o.m,'gravity',o.gravity);
            %[o.KK0,o.KKg] = fGeneralizedStiffnessMatrixBeam(o.s_span,o.EI,o.PhiK,o.bOrth,o.PhiV,o.m,o.Mtop,o.gravity,[],[],bStiffeningMtop,bStiffeningSelfWeight)

            if o.bStiffnessFromGM
                if isempty(o.ModeOmegas); error('ModeOmegas needs to be set when using bStiffnessFromGM'); end
                GM    = diag( o.MM(7:end,7:end));
                GK    = o.ModeOmegas(:).^2 .* GM(:);
                KKf = diag(GK);
                o.KK = zeros(6+o.nf,6+o.nf);
                o.KK(7:end,7:end)=KKf;
                o.KK0=o.KK;
            else
                o.KK=o.KK0;
                %o.KK=o.KK0+o.KKg; % TODO
            end
        end

    end

    function o=computeDampingMatrix(o)
        if isequal(o.DampingMethod,'Rayleigh')
            Rayleigh_beta  = o.DampingParams(1);
            Rayleigh_alpha = o.DampingParams(2);
            o.DD=Rayleigh_beta*o.KK+Rayleigh_alpha*o.MM;
        elseif isequal(o.DampingMethod,'StiffnessProportional')
            disp('[WARN] StiffnessProportional only for Orthogonal Modes');
            o.DD = zeros(6+o.nf,6+o.nf);
            DDf=zeros(o.nf,o.nf);
            for j=1:o.nf
                LogDec = o.DampingParams(1);
                Omega=sqrt(o.KK(6+j,6+j)./o.MM(6+j,6+j));
                DDf(j,j)=LogDec./(pi*Omega).*o.KK(6+j,6+j) ;
            end
            o.DD(7:end,7:end)= DDf;
        elseif isequal(o.DampingMethod,'None')
            o.DD = zeros(6+o.nf,6+o.nf);
        else
            error('Unknown damping method %s',o.DampingMethod);
        end
    end

    % --------------------------------------------------------------------------------}
    %% ---  Connection
    % --------------------------------------------------------------------------------{
    function connectTo(o,Child,sPoint,MyBodyPoint,varargin)

        if isequal(o.Type,'Rigid') || isequal(o.Type,'Ground')
            if ischar(MyBodyPoint); error('Connection Point cannot be character for rigid body'); end;
            s_C_inB=MyBodyPoint;
            c=Connection('Name',[o.Name '-' Child.Name],'RelPoint',s_C_inB,varargin{:});
        else
            if ischar(MyBodyPoint) 
                switch MyBodyPoint
                    case 'FirstPoint'
                        i_C_inB=1;
                    case 'LastPoint'
                        i_C_inB=o.nSpan;
                    otherwise
                        error('Unsupported Body point option %s',MyBodyPoint)
                end
            else
                % trying to find index
                error('todo')
            end
            s_C_inB=o.s_P(:,i_C_inB);
            c=Connection('Name',[o.Name '-' Child.Name],'RelPoint',s_C_inB,'ParentNode',i_C_inB,varargin{:});
        end
        if isempty(o.Children   ); o.Children    = cell(0); end;
        if isempty(o.Connections); o.Connections = cell(0); end;
        for i=1:length(o.Children)
            % octave comment
            %if o.Children{i} == Child % check for handles
            %    error('Cannot connect twice to the same body');
            %end
        end
        o.Children{end+1}=Child;
        o.Connections{end+1}=c;
    end

    function n=setupDOFIndex(o,n)
        nForMe=o.nf;
        % Setting my dof index
        o.I_DOF=n+ (1:nForMe); 
        % Update
        n=n+nForMe;
        for ic=1:length(o.Connections)
            % --- Connection first
            nForConn=o.Connections{ic}.nj;
            o.Connections{ic}.I_DOF=n+(1:nForConn);
            % Update
            n=n+nForConn;
            % --- Then Children
            n=o.Children{ic}.setupDOFIndex(n);
        end
    end

    function updateChildrenKinematicsNonRecursive(p,q,qdot,qddot)
        % At this stage all the kinematics of the body p are known


        % Useful variables
        R_0p =  p.R_0b;
        B_p  =  p.B;
        r_0p  = p.r_O;

        nf_all_children=0; % TODO not generic, we need to look at I_DOF
        for ic=1:length(p.Children)
           nf_all_children=nf_all_children+p.Children{ic}.nf;
        end

        for ic=1:length(p.Children)
            i    = p.Children{ic}   ;
            conn = p.Connections{ic};
            % Position of connection point in P and 0 system
            r_pi_inP= conn.s_C_inB;
            r_pi    = R_0p * r_pi_inP; 
            % Flexible influence to connection point
            if isequal(p.Type,'FlexibleBeam')
                % TODO make this general...
                if conn.ParentNode==p.nSpan
                    iNode=conn.ParentNode;
                    %CyT= - [Twr.PhiV{1}(3,end) Twr.PhiV{2}(3,end)];
                    alpha_y= - p.V(3,iNode);
                    alpha_z=   p.V(2,iNode);
                    R_pc=fRoty(alpha_y)*fRotz(alpha_z);


                    Bx_pc=zeros(3,p.nf);
                    Bt_pc=zeros(3,p.nf);
                    for j=1:p.nf
                        Bx_pc(:,j)=p.PhiU{j}(:,iNode);
                        Bt_pc(:,j)=[0; -p.PhiV{j}(3,iNode); p.PhiV{j}(2,iNode)];
                    end
                else
                    error('B matrix with other connection point')
                end
            else
                R_pc=eye(3);
                Bx_pc=zeros(3,0);
                Bt_pc=zeros(3,0);
            end
            % Joint influence to next body (R_ci, B_ci)
            conn.updateKinematics(q,qdot);
%             if ic==2
%             end

            % Full connection p and j
            R_pi   = R_pc*conn.R_ci  ;
            Bx_pi  = [Bx_pc R_pc*conn.B_ci(1:3,:)];
            Bt_pi  = [Bt_pc R_pc*conn.B_ci(4:6,:)];


            % Rotation of body i is rotation due to p and j
            R_0i = R_0p * R_pi;

            B_i=fBMatRecursion(B_p,[Bx_pi; Bt_pi],R_0p,r_pi);
            B_i_inI = [R_0i' * B_i(1:3,:);  R_0i' * B_i(4:6,:)];
            nf_i=i.nf;
            % TODO  not generic, we need to look at I_DOF, and this assumes all child has same amount of dofs...
            I_B = zeros(nf_i,nf_all_children);
            iOff=(ic-1)*nf_i;
            I_B(1:nf_i,iOff+(1:nf_i))=eye(nf_i);
            BB_i_inI=[B_i_inI zeros(6,nf_all_children) ; zeros(nf_i, size(B_i_inI,2)) I_B];

%             BB_i_inI=[B_i_inI zeros(6,nf_i)            ; zeros(nf_i, size(B_i_inI,2)) eye(nf_i)];
            i.B     = B_i    ;
            i.B_inB = B_i_inI;
            i.BB_inB = BB_i_inI;

            % Index of relevant DOFs
            n_to_i    = size(B_i,2)         ; % Better be correct..
            IHat      = [p.I_DOF conn.I_DOF];
            I_to_i    = 1:n_to_i            ;
            I_with_i  = 1:(n_to_i+nf_all_children); % TODO need account of I_DOF

            % --- Updating Position of child body 
            r_0i = r_0p + r_pi  ; % in 0 system
            gz   = q    (i.I_DOF);
            gzp  = qdot (i.I_DOF);
            gzpp = qddot(i.I_DOF);

            % Velocity
            v_pi_inP  = Bx_pi*qdot(IHat)    ;
            om_pi_inP = Bt_pi*qdot(IHat)    ;

            % All velocities
            %v_0(1:3)      = B_i(1:3,:) * qdot(1:n_to_i);
            %v_0(4:6)      = B_i(4:6,:) * qdot(1:n_to_i);
            %v_0(7:6+i.nf) = qdot(n_to_i+1:i.nf);
            v_i_in0 = BB_i_inI * qdot(I_with_i);
            %  
            % V-Accelerations in P and 0
            a_i_v_in0=zeros(6+i.nf,1); % All v accelerations
            a_i_v_inP =p.a_O_v_inB + ...
                     cross(p.om_O_inB,cross(p.om_O_inB, r_pi_inP) )+...
                     2*cross(p.om_O_inB,  v_pi_inP ) + ...
                     cross(p.omp_O_v_inB, r_pi_inP);

            omp_i_v_inP =p.omp_O_v_inB +  cross(p.om_O_inB, om_pi_inP);
            omp_i_v_inI_bis = R_pi'*omp_i_v_inP; % <<<<<<<< ALTERNATIVE
            omp_i_v_inI =zeros(3,1);
            dom=zeros(3,1);
            om=zeros(3,1);
            for j=size(B_i_inI,2):-1:1
                dom         = B_i_inI(4:6,j)*qdot(j)      ;
                omp_i_v_inI = omp_i_v_inI +  cross(dom,om);
                om          = om + dom                    ;
            end
            om_in0     = B_i(4:6,:) * qdot(1:n_to_i);
            om_inI_bis = R_0i *om_in0; % <<<<<<<< ALTERNATIVE

            % Accelerations in 0
            a_i_v_in0(1:3)=R_0p*a_i_v_inP;
            a_i_v_in0(4:6)=R_0p*omp_i_v_inP;
            a_i_v_in0(4:6)=R_0i*omp_i_v_inI; % <<<<<<<< ALTERNTIVE
            a_i_v_in0(7:6+i.nf)=gzpp;


            i.R_pb = R_pi ;
            i.updateKinematics(r_0i,R_0i,gz,v_i_in0,a_i_v_in0)
            if ~isequal(i.Type,'Rigid')
                if i.nf>0
                    i.computeMassMatrix(); % Triggers 
                    i.computeInertiaForces();
                end
            end
%             nf_N=0;
%             BB_N_inN =[B_N_inN zeros(6,nf_N) ; zeros(nf_N, size(B_N_inN,2)) eye(nf_N)];
%             % Update of mass matrix
%             MM_N= BB_N_inN'*Nac.MM*BB_N_inN;
%             KK_N= BB_N_inN'*Nac.KK*BB_N_inN;



        end
    end

    function MqB=getFullMContrib(o)
        MqB=o.BB_inB'*o.MM*o.BB_inB;
    end

    function M=getFullM(o,M)
        if ~isequal(o.Type,'Ground')
            MqB=o.getFullMContrib();
            n=size(MqB,1);
            M(1:n,1:n)=M(1:n,1:n)+MqB;
        end
        for ic=1:length(o.Children)
             M=o.Children{ic}.getFullM(M);
         end
    end
    function K=getFullK(o,K)
        if ~isequal(o.Type,'Ground')
            KqB=o.BB_inB'*o.KK*o.BB_inB;
            n=size(KqB,1);
            K(1:n,1:n)=K(1:n,1:n)+KqB;
        end
        for ic=1:length(o.Children)
             K=o.Children{ic}.getFullK(K);
         end
    end



    % --------------------------------------------------------------------------------}
    %% ---  
    % --------------------------------------------------------------------------------{
    function checkModesOrthogonality(o)
        nModes=length(o.PhiU);
        ModeOrth=zeros(nModes,nModes);
        if nModes>0
            for j=1:nModes
                for i=j:nModes
                    ModeOrth(i,j)=trapz(o.s_span,o.PhiU{i}(1,:).*o.PhiU{j}(1,:)+ o.PhiU{i}(2,:).*o.PhiU{j}(2,:)+ o.PhiU{i}(3,:).*o.PhiU{j}(3,:));
                end
            end
            Cross=triu(ModeOrth,1);
        end
        %fprintf('ModeOrth: max %.2e   avg %.2e\n',max(Cross(:)),mean(Cross(:)));
    end




    function v=getMeanLine(o)
        r_C = o.b2g_s(o.s_C_inB);
        % Vector 
        v=[o.r_O(:) r_C(:)];
    end

    function plotModes(o)
        for j=1:length(o.PhiU)
            figure
            title(sprintf('Mode %d',j));
            plot(o.s_span, o.PhiU{j});
            legend('x','y','z')
            hold all
        end
    end

    function r=b2g_s(o,s)
        r = o.r_O + o.R_0b*s;
    end

    % --------------------------------------------------------------------------------}
    %% --- Flex 
    % --------------------------------------------------------------------------------{
    function setCompatible(o,b)
        o.bCompatibility=b;
        if b
            o.bOrthogonal=b;
            o.bStiffnessFromGM=b;
        end
    end
    % --------------------------------------------------------------------------------}
    %% --- Helper functions 
    % --------------------------------------------------------------------------------{
    function check_not_empty(o,varargin)
        for i=1:length(varargin)
            if isempty(o.(varargin{i}));
                error('Property %s needs to be set',varargin{i});
            end
        end
    end
    function check_same_size(o,varargin)
        S=size(o.(varargin{1}));
        for i=2:length(varargin)
            if ~isequal(S, size(o.(varargin{i})))
                error('Property %s does not have the same dimension as %s',varargin{i},varargin{1});
            end
        end
    end
    function line_vector_me(o,varargin)
        for i=1:length(varargin)
            o.(varargin{i})=o.(varargin{i})(:)';
        end
    end


end % methods

end % class
