function [MM] = fGMBeamStandaloneSimplified(s_G,s_span,m,U)
% Performing full integration of mass matrix without shape integral functions
% for a beam
% 
% INPUTS
%  - s_G    : [m] 3 x nSpan, location of the cross sections COG points, when the beam is undeflected
%  - s_span : [m] span integration variable (e.g. s_G(1,:) if beam along x, or s_G(3,:) if beam along z)
%  - m      : [kg/m] cross section mass along the beam
%
%
% OPTIONAL INPUTS:
%  - U: cell of shape functions (displacements functions)
%       Each shape function U{j}, is a 3 x nSpan matrix corresponding 
%         to the displacements of the center of mass s_G
%       If omitted, then rigid body (6x6) mass matrix is returned
%        
%
% AUTHOR: E. Branlard
%
% --- Optional arguments
if ~exist('U'         ,'var'); U          = []      ; end
ng =length(U);

% --- local variables for fast integration
nSpan = length(s_span);
% dxx = diff(s_span)';
IW  =zeros(1,nSpan); 
for i=1:nSpan-1 
    L           = s_span(i+1) - s_span(i);
    IW  (i)     = IW(i) + L/2          ;
    IW  (i+1)   = L/2                        ;
end

% Mass integration
Mass = trapzs(s_span, m);

% --- Mxt = -\int [~s] dm    =  -Skew(sigma+Psi g)
% S(1:3)=trapz(s_span, s_G(1:3,:) .* dm  , 2) ;
% Mxt = -fSkew(S);
S(1)=trapzs(s_span, s_G(1,:) .* m);
S(2)=trapzs(s_span, s_G(2,:) .* m);
S(3)=trapzs(s_span, s_G(3,:) .* m);
Mxt=[0    S(3)  -S(2) ;
    -S(3)  0     S(1) ; 
    +S(2) -S(1)   0  ];
% --- Mtt = - \int [~s][~s] dm
s11=trapzs(s_span, s_G(1,:) .* s_G(1,:) .* m);
s12=trapzs(s_span, s_G(1,:) .* s_G(2,:) .* m);
s13=trapzs(s_span, s_G(1,:) .* s_G(3,:) .* m);
s22=trapzs(s_span, s_G(2,:) .* s_G(2,:) .* m);
s23=trapzs(s_span, s_G(2,:) .* s_G(3,:) .* m);
s33=trapzs(s_span, s_G(3,:) .* s_G(3,:) .* m);
Mtt=zeros(3,3);
Mtt(1,1)= s22+s33; Mtt(1,2)= -s12   ; Mtt(1,3)= -s13   ;
Mtt(2,1)= -s12   ; Mtt(2,2)= s11+s33; Mtt(2,3)= -s23   ;
Mtt(3,1)= -s13   ; Mtt(3,2)= -s23   ; Mtt(3,3)= s11+s22;

% --- Mxg = \int Phi dm  =   Psi
Mxg=NaN(3,ng);
for j=1:ng
    Mxg(1,j) = trapzs(s_span, U{j}(1,:) .* m);
    Mxg(2,j) = trapzs(s_span, U{j}(2,:) .* m);
    Mxg(3,j) = trapzs(s_span, U{j}(3,:) .* m);
end

% --- Mtg  = \int [~s] Phi dm  
% With torsion contribution if any
Mtg=zeros(3,ng);
for j=1:ng
    Mtg(1,j)= trapzs(s_span,  (- s_G(3,:) .* U{j}(2,:) + s_G(2,:) .* U{j}(3,:)).*m );
    Mtg(2,j)= trapzs(s_span,  (+ s_G(3,:) .* U{j}(1,:) - s_G(1,:) .* U{j}(3,:)).*m );
    Mtg(3,j)= trapzs(s_span,  (- s_G(2,:) .* U{j}(1,:) + s_G(1,:) .* U{j}(2,:)).*m );
end

% --- Mgg  = \int Phi^t Phi dm  =  Sum Upsilon_kl(i,i)
Mgg=zeros(ng,ng);
for i=1:ng; for j=1:ng
    Mgg(i,j)=trapzs(s_span, (U{i}(1,:).*U{j}(1,:) + U{i}(2,:).*U{j}(2,:) + U{i}(3,:).*U{j}(3,:) ).*m);
end; end;

% --- Building M and making it symmetric
MM=zeros(6+ng, 6+ng);
MM(1:3      , 1:3)      = Mass*eye(3); % Mxx
MM(1:3      , 4:6)      = Mxt;
MM(1:3      , 6+(1:ng)) = Mxg;
MM(4:6      , 4:6)      = Mtt;
MM(4:6      , 6+(1:ng)) = Mtg;
MM(6+(1:ng) , 6+(1:ng)) = Mgg;
MM=triu(MM)+triu(MM,1)';



function zz = trapzs(~,yy)
    % yy needs to be a line vector
    %zz =  (yy(1:nSpan-1) + yy(2:nSpan))/2 * dxx;
    zz = sum(yy.*IW);
end % trapz1d
end
