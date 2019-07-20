function [MM,MMG_d0]=fGMBeamCompatilibity(gzf,s_G0,s_P,rho_G0,s_span,m,jxxG,PhiUG,PhiV,V_tot, Method,V0)
%% Documentation
%   Computes  Mass Matrix with Flex-like compatibility
%   For that you need (correctly): 
%     - Axial correction
%     - Integration weight IW_xm
%     - UG shape functions for the "mode" part of the matrix
%   For that you also need (Unfortunately): 
%     - Some of the terms computed at G_d0 and not G 
% 
% NOTE: 
% -  Methods: ShapeIntegrals_UG_G0 and ShapeIntegralsStandalone_UG_G0 are IDENTICAL!!!
%
%  Differences between points:
%    s_G_d0 = s_P + rho_0 =  s_P0 + gj PhiU  + rho_0
%    s_G0   =                s_P0 +          + rho_0
%    s_G_lin=                s_P0 + gj PhiUG + rho_0
%    s_G    =                s_P0 + gj PhiUG + rho


bOrth=true;

% --- Span and Mass Integration weights  IW and IW_x 
nSpan=length(s_span);
IW =zeros(1,nSpan); IW_x =zeros(1,nSpan); IW_xm=zeros(1,nSpan); Imom = 0;
for I=1:nSpan-1 
    L         = s_span(I+1) - s_span(I);
    IW  (I)     = IW(I) + L/2          ;
    IW  (I+1)   = L/2                        ;
    IW_x(I)   = IW_x(I) + (s_span(I)/2 + L/6)*L;
    IW_x(I+1) = (s_span(I)/2 + L/3)*L;
end
for I=1:nSpan
    IW_xm(I) = IW_x(I)*m(I);
    Imom  = Imom  + IW_xm(I)*s_span(I); % Imom = \int x^2.m(x) dx
end


% s_G_d0 = s_P0 + gj PhiU  + rho_0
s_G_d0 = s_P;
s_G_d0(2,:) = s_G_d0(2,:)+rho_G0(2,:); 
s_G_d0(3,:) = s_G_d0(3,:)+rho_G0(3,:); % TODO torsion effects


%% 
if isequal(Method,'ShapeIntegrals_UG_G0')
    if ~exist('V0','var'); error('Please provide V0 with this method'); end;
    [Psi,Upsilon_kl,Sigma_kl,sigma_kl,sigma,Mass,ImomX,I_Jxx,GMJxx,M1] = fShapeIntegralsBeam(s_G0,s_span,m,jxxG,PhiUG,PhiV,'IW_xm',IW_xm,'bOrth',bOrth,'bAxialCorr',true,'V0',V0);
    [MMG_SI] = fGMBeamShapeIntegrals(gzf,Psi,Upsilon_kl,Sigma_kl,sigma_kl,sigma,Mass,ImomX,I_Jxx,GMJxx,bOrth,true,M1);

elseif isequal(Method,'ShapeIntegralsStandalone_UG_G0')
    if ~exist('V0','var'); error('Please provide V0 with this method'); end;
    [MMG_SI  ]= fGMBeamStandaloneShapeIntegrals(s_G0,s_span,m,jxxG,PhiUG,PhiV,gzf,'IW_xm',IW_xm,'bOrth',bOrth,'bAxialCorr',true,'V0',V0);

elseif isequal(Method,'DirectIntegration_UG_G_d0')
    % MMG_d0, but required no matter what...
else
    error('Unknown method %s',Method)
end

[MMG_d0]= fGMBeamStandalone              (s_G_d0 ,s_span,m,jxxG,PhiUG,PhiV,IW,IW_xm,bOrth,true,V_tot);

%% HACK to get results like Flex
if isequal(Method,'ShapeIntegrals_UG_G0') || isequal(Method,'ShapeIntegralsStandalone_UG_G0')
    MMH=MMG_SI;
    nDOF=size(MMH,2);
    % Mtf compatibility
    I = 4:6 ; J = 7:nDOF;
    MMH(4,J)=MMG_d0(4,J);
    MMH(J,I)=MMH(I,J)';
    % Mxt compatiblity - Only first line changed
    MMH(2,4)=MMG_d0(2,4); MMH(4,2)=MMG_d0(2,4);
    MMH(3,4)=MMG_d0(3,4); MMH(4,3)=MMG_d0(3,4); % important when you have some BetaC!
    % Mtt compatbility - 4 terms and constant 5-6
    MMH(4,4)=MMG_d0(4,4);
    MMH(5,4)=MMG_d0(5,4); MMH(4,5)=MMH(5,4);
    MMH(6,4)=MMG_d0(6,4); MMH(4,6)=MMH(6,4);
    MMH(5,6)=0; MMH(6,5)=0;
    MMH(6,6)=Imom;
    MMH(5,5)=MMH(6,6);
    % We return
    MM=MMH;
elseif isequal(Method,'DirectIntegration_UG_G_d0')
    MMH_d0=MMG_d0;
    % HACK compatibility - constant 5-6
    MMH_d0(5,6)=0; MMH_d0(6,5)=0;
    MMH_d0(6,6)=Imom; % Annoying
    MMH_d0(5,5)=MMH_d0(6,6);
    % We return
    MM=MMH_d0;
end




