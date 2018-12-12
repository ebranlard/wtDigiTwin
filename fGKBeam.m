function [KK0,KKg] = fGKBeam(s_span,GKt,EIy,EIz,K,V,V0,K0,varargin);
%     ,V,m,Mtop,gravity,xLumped,MLumped,bStiffeningMtop,bStiffeningSelfWeight)
% 
% OUTPUTS
%   KK0: Stiffness matrix without geometrical stiffening
%   KKg: Geometrical stiffeness matrix. 
%   The total stiffness matrix should then be KK0+KKg
%
% AUTHOR: E.Branlard

% --- Optional arguments
%p=inputParser();
if mod(length(varargin),2)~=0; error('Use pairs of key/value for optional arguments'); end
p=fInputParser();
%p.KeepUnmatched=true;
p.addParameter('bOrth',false);
% For stiffning
p.addParameter('Mtop'                 ,0       ,@isnumeric);
p.addParameter('m'                    ,s_span*0,@isnumeric);
p.addParameter('gravity'              ,9.81    ,@isnumeric);
p.addParameter('bStiffeningMtop'      ,true    )           ;
p.addParameter('bStiffeningSelfWeight',true    )           ;
p.parse(varargin{:}); p=p.Results;

% --- Optional arguments for Stiffening correction
% if ~exist('xLumped'              ,'var'); xLumped               = []  ; end;
% if ~exist('MLumped'              ,'var'); MLumped               = []  ; end;
% 

% --- Derived parameters
nSpan=length(s_span);
nf=length(K); % number of modes
bStiffening=(p.bStiffeningMtop & ~isempty(p.Mtop)) | (p.bStiffeningSelfWeight & ~isempty(p.m));
% 
% if ~isempty(xLumped) ||  ~isempty(MLumped)
%     error('So far this script is only for Mtop, Mlumped to be implemented')
% end

% --- Init
KK0 = zeros(6+nf,6+nf);
KKg = zeros(6+nf,6+nf);

% --------------------------------------------------------------------------------}
%% --- Stiffening  
% --------------------------------------------------------------------------------{
if bStiffening
    %% --- Axial force 
    Pacc_SW = fcumtrapzlr(s_span, -p.m * p.gravity) ;
    Pacc_MT = -p.Mtop * p.gravity*ones(1,nSpan);
    Pacc    = zeros(1,nSpan)               ;
    % TopMass contribution to Pacc
    if p.bStiffeningMtop
        Pacc=Pacc+Pacc_MT;
    end
    if p.bStiffeningSelfWeight
        Pacc=Pacc+Pacc_SW;
    end

    % Method 2
    KKCorr = zeros(nf,nf)  ;
    for i=1:nf
        for j=1:nf
            % TODO maybe V0 should be put there, but I don't think so?
            %xx=trapz(s_span, Pacc .* V{i}(1,:).* o.V{j}(1,:));
            yy=trapz(s_span, Pacc .*V{i}(2,:) .* V{j}(2,:)); 
            zz=trapz(s_span, Pacc .*V{i}(3,:) .* V{j}(3,:));
            KKCorr(i,j)=yy+zz;
        end
    end
    KKg(7:end,7:end)= KKCorr;
end


% --- Transformation matrix from blade to section
R_sb=cell(1,nSpan);
for i=1:length(s_span)
%     theta_x=-Bld.BetaC(i); theta_y=-V0(3,i); theta_z= V0(2,i);
    theta_x= V0(1,i);  % Pincipal axis "twist"
    theta_y=-V0(3,i);  % duz/dx = -\theta_y
    theta_z= V0(2,i);  % duy/dx =  \theta_z
%     R_sb{i}=fTmatrix(2,3,1,cos(theta_y), sin(theta_y),cos(theta_z),sin(theta_z),cos(theta_x), sin(theta_x));
    R_sb{i}=fRotx(theta_x)'*fRotz(theta_z)'* fRoty(theta_y)';
end


KKf=zeros(nf,nf);
% --- Integrating to get generalized stiffness
for i=1:nf
    K_si=zeros(3,nSpan);
    K_bi=zeros(3,nSpan);
    % Mode curvature plus correction due to torsion/prebend
    K_bi(1,:) =   K{i}(1,:);
    K_bi(2,:) =   K{i}(2,:) + (V0(3,:).*K{i}(1,:)+V{i}(1,:).*K0(3,:) ) ;
    K_bi(3,:) =   K{i}(3,:) - (V0(2,:).*K{i}(1,:)+V{i}(1,:).*K0(2,:) ) ; 
    for j=1:nf
        K_sj=zeros(3,nSpan);
        K_bj=zeros(3,nSpan);
        % Full Curvature = Mode curvature + Correction due to torsion/prebend interaction
        K_bj(1,:) =   K{j}(1,:);
        K_bj(2,:) =   K{j}(2,:) + (V0(3,:).*K{j}(1,:)+V{j}(1,:).*K0(3,:) ) ;
        K_bj(3,:) =   K{j}(3,:) - (V0(2,:).*K{j}(1,:)+V{j}(1,:).*K0(2,:) ) ; 
        % Transforming curvature from blade to section coordinate
        %     K_si=K_bi;
        %     K_sj=K_bj;
        for is=1:length(s_span)
            K_si(:,is)=R_sb{is}*K_bi(:,is);
            K_sj(:,is)=R_sb{is}*K_bj(:,is);
        end

        % Integration EI*ki*kj
        %xx=trapz(s_span, EIx .* K{i}(1,:).* K{j}(1,:));
        %yy=trapz(s_span, EIy .* K{i}(2,:).* K{j}(2,:));
        %zz=trapz(s_span, EIz .* K{i}(3,:).* K{j}(3,:));
        xx = trapz(s_span, GKt.*K_si(1,:).*K_sj(1,:) );
        yy = trapz(s_span, EIy.*K_si(3,:).*K_sj(3,:) ); % NOTE:  My = -EIy k_z
        zz = trapz(s_span, EIz.*K_si(2,:).*K_sj(2,:) ); % NOTE:  Mz =  EIz k_y
        %
        KKf(i,j)=xx+yy+zz;
    end
end

if p.bOrth
    KKf=KKf.*eye(nf); % Ensuring orthogonality
end
KK0(7:end,7:end)= KKf;
