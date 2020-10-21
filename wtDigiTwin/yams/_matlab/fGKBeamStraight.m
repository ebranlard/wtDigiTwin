function [KK0,KKg] = fGKBeamStraight(s_span,EIx,EIy,EIz,PhiK,bOrth,PhiV,m,Mtop,gravity,xLumped,MLumped,bStiffeningMtop,bStiffeningSelfWeight)
% 
% OUTPUTS
%   KK0: Stiffness matrix without geometrical stiffening
%   KKg: Geometrical stiffeness matrix. 
%   The total stiffness matrix should then be KK0+KKg

% --- Optional arguments for Stiffening correction
if ~exist('m'                    ,'var'); m                     = []  ; end
if ~exist('PhiV'                 ,'var'); PhiV                  = []  ; end
if ~exist('gravity'              ,'var'); gravity               = 9.81; end;
if ~exist('Mtop'                 ,'var'); Mtop                  = []  ; end;
if ~exist('xLumped'              ,'var'); xLumped               = []  ; end;
if ~exist('MLumped'              ,'var'); MLumped               = []  ; end;
if ~exist('bStiffeningMtop'      ,'var'); bStiffeningMtop       = true; end
if ~exist('bStiffeningSelfWeight','var'); bStiffeningSelfWeight = true; end;


% --- Derived parameters
nf=length(PhiK); % number of modes
bStiffening=(bStiffeningMtop & ~isempty(Mtop)) | (bStiffeningSelfWeight & ~isempty(m));

if ~isempty(xLumped) ||  ~isempty(MLumped)
    error('So far this script is only for Mtop, Mlumped to be implemented')
end

% --- Init
KK0 = zeros(6+nf,6+nf);
KKg = zeros(6+nf,6+nf);

% --------------------------------------------------------------------------------}
%% --- Stiffening  
% --------------------------------------------------------------------------------{
if bStiffening
    nSpan=length(s_span);
    %% --- Axial force 
%     keyboad
    Pacc_SW = fcumtrapzlr(s_span, -m * gravity) ;
    Pacc_MT = -Mtop * gravity*ones(1,nSpan);
    Pacc    = zeros(1,nSpan)               ;
    % TopMass contribution to Pacc
    if bStiffeningMtop
        Pacc=Pacc+Pacc_MT;
    end
    if bStiffeningSelfWeight
        Pacc=Pacc+Pacc_SW;
    end

    % Method 2
    KKCorr = zeros(nf,nf)  ;
    for i=1:nf
        for j=1:nf
            %xx=trapz(s_span, Pacc .* PhiV{i}(1,:).* o.PhiV{j}(1,:));
            yy=trapz(s_span, Pacc .*PhiV{i}(2,:) .* PhiV{j}(2,:));
            zz=trapz(s_span, Pacc .*PhiV{i}(3,:) .* PhiV{j}(3,:));
            KKCorr(i,j)=yy+zz;
        end
    end
    KKg(7:end,7:end)= KKCorr;
end



%% --- Stiffness matrix 
KKf = zeros(nf,nf)  ;
for i=1:nf
    for j=1:nf
        xx=trapz(s_span, EIx .* PhiK{i}(1,:).* PhiK{j}(1,:)); % TODO TODO
        yy=trapz(s_span, EIy .* PhiK{i}(2,:).* PhiK{j}(2,:));
        zz=trapz(s_span, EIz .* PhiK{i}(3,:).* PhiK{j}(3,:));
        KKf(i,j)=xx+yy+zz;
    end
end
if bOrth
    KKf=KKf.*eye(nf); % Ensuring orthogonality
end
KK0(7:end,7:end)= KKf;
