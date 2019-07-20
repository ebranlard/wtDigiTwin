function [MM,Jxx,I_Jxx,GMJxx] = fGMBeamStandalone(s_G,s_span,m,jxxG,U,V,IW,IW_xm,bOrth,bAxialCorr,V_tot,Peq_tot,IW_U)
% Performing full integration of mass matrix without shape integral functions
% NOTE: Beam assumed to be along x for now (only because of Jxx)
% 
% INPUTS
%  - s_G    : [m] 3 x nSpan , location of cross sections COG
%  - s_span : [m] span integration variable (e.g. s_G(1,:))
%  - m      : [kg/m] cross section mass along the beam
%  - jxxG   : [kg.m] second moment of inertia of cross section
%
%
% OPTIONAL INPUTS:
%  - JxxG, if omitted, assumed to be 0
%  - U , if omitted, then rigid body (6x6) mass matrix is returned
%
% AUTHOR: E. Branlard
%
% --- Optional arguments
if ~exist('jxxG'      ,'var'); jxxG       = 0       ; end;
if ~exist('U'         ,'var'); U          = []      ; end
if ~exist('V'         ,'var'); V          = []      ; end
if ~exist('bOrth'     ,'var'); bOrth      = false   ; end
if ~exist('bAxialCorr','var'); bAxialCorr = false   ; end
if ~exist('IW_U','var')      ; IW_U             = []; end
if ~exist('Peq_tot','var')   ; Peq_tot          = []; end
if ~exist('V_tot'  ,'var')   ; V_tot          = []; end
nf =length(U);

% --- Safety checks
if ~isempty(V)
    if length(V)~=nf; error('Size of V should match size of U'); end;
end

% --- local variables for fast integration
nSpan = length(s_span);
% dxx = diff(s_span)';
if isempty(IW) 
    IW  =zeros(1,nSpan); 
    for i=1:nSpan-1 
        L           = s_span(i+1) - s_span(i);
        IW  (i)     = IW(i) + L/2          ;
        IW  (i+1)   = L/2                        ;
    end
end




% Mass integration
Mass = trapzs(s_span, m);

% --- Torsion-related variables - Zeros by default
Jxx = trapzs(s_span, jxxG); % Imomx
GMJxx=zeros(1,nf);
I_Jxx=zeros(1,nf);
for j=1:length(V) % NOTE: only computed if V not empty!
    VJ=jxxG.*V{j}(1,:);
    GMJxx(j)=trapzs(s_span,V{j}(1,:).*VJ);
    %I_Jxx(j)=sum(Bld.GFTJ.*Bld.VEVx(j,:)); 
    I_Jxx(j)=trapzs(s_span,VJ);
end


% --- Mxt = -\int [~s] dm    =  -Skew(sigma+Psi g)
% S(1:3)=trapz(s_span, s_G(1:3,:) .* dm  , 2) ;
% Mxt = -fSkew(S);
S(1)=trapzs(s_span, s_G(1,:) .* m);
S(2)=trapzs(s_span, s_G(2,:) .* m);
S(3)=trapzs(s_span, s_G(3,:) .* m);
Mxt=[0    S(3)  -S(2) ;
    -S(3)  0     S(1) ; 
    +S(2) -S(1)   0  ];

if bAxialCorr
    % --- Variables for axial correction
    % FT=fcumtrapzlr(s_span,m);
    FT = - fliplr( cumtrapz( fliplr(s_span), fliplr(m)) ) ;

    if isempty(V_tot); error('Please provide Vtot for axial correction'); end
    m15=+trapzs( s_span, V_tot(3,:).*FT);
    m16=-trapzs( s_span, V_tot(2,:).*FT);
    Mxt(1,2)=m15;
    Mxt(1,3)=m16;
end

% --- Mtt = - \int [~s][~s] dm
if ~isempty(IW_xm) 
    s11= sum(          IW_xm .* s_G(1,:)     );
    s12= sum(          IW_xm .* s_G(2,:)     );
    s13= sum(          IW_xm .* s_G(3,:)     );
else
    s11=trapzs(s_span, s_G(1,:) .* s_G(1,:) .* m);
    s12=trapzs(s_span, s_G(1,:) .* s_G(2,:) .* m);
    s13=trapzs(s_span, s_G(1,:) .* s_G(3,:) .* m);
end
s22=trapzs(s_span, s_G(2,:) .* s_G(2,:) .* m);
s23=trapzs(s_span, s_G(2,:) .* s_G(3,:) .* m);
s33=trapzs(s_span, s_G(3,:) .* s_G(3,:) .* m);
Mtt=zeros(3,3);
Mtt(1,1)= s22+s33; Mtt(1,2)= -s12   ; Mtt(1,3)= -s13   ;
Mtt(2,1)= -s12   ; Mtt(2,2)= s11+s33; Mtt(2,3)= -s23   ;
Mtt(3,1)= -s13   ; Mtt(3,2)= -s23   ; Mtt(3,3)= s11+s22;
% Torsion contribution
Mtt(1,1)=Mtt(1,1)+Jxx; 

% --- Mxg = \int Phi dm  =   Psi
Mxg=NaN(3,nf);
for j=1:nf
    Mxg(1,j) = trapzs(s_span, U{j}(1,:) .* m);
    Mxg(2,j) = trapzs(s_span, U{j}(2,:) .* m);
    Mxg(3,j) = trapzs(s_span, U{j}(3,:) .* m);
end
if bAxialCorr
    if ~isempty(V_tot) && ~isempty(Peq_tot); error('Provide either V_tot or Peq_tot'); end;
    if ~isempty(V_tot)
        for j=1:nf
            Mxg(1,j)= trapzs(s_span, -V{j}(2,:).*V_tot(2,:).*FT  - V{j}(3,:).*V_tot(3,:).*FT); 
        end
    elseif ~isempty(Peq_tot)
        for j=1:nf
            Mxg(1,j) = trapzs(s_span, U{j}(2,:).*Peq_tot(2,:) + U{j}(3,:).*Peq_tot(3,:) );
        end
    else
        error('Please provide Vtot of Peq_tot for axial correction');
    end
end



% --- Mtg  = \int [~s] Phi dm  
% With torsion contribution if any
Mtg=zeros(3,nf);
if ~isempty(IW_xm) 
    for j=1:nf
        Mtg(1,j)= trapzs(s_span, (- s_G(3,:) .* U{j}(2,:)  + s_G(2,:) .* U{j}(3,:)).*m  ) + I_Jxx(j);
        Mtg(2,j)= trapzs(s_span, ( + s_G(3,:) .* U{j}(1,:).*m))  - sum(IW_xm.* U{j}(3,:));
        Mtg(3,j)= trapzs(s_span, ( - s_G(2,:) .* U{j}(1,:).*m))  + sum(IW_xm.* U{j}(2,:));
    end
else
    for j=1:nf
        Mtg(1,j)= trapzs(s_span,  (- s_G(3,:) .* U{j}(2,:) + s_G(2,:) .* U{j}(3,:)).*m ) + I_Jxx(j);
        Mtg(2,j)= trapzs(s_span,  (+ s_G(3,:) .* U{j}(1,:) - s_G(1,:) .* U{j}(3,:)).*m );
        Mtg(3,j)= trapzs(s_span,  (- s_G(2,:) .* U{j}(1,:) + s_G(1,:) .* U{j}(2,:)).*m );
    end
end

% --- Mgg  = \int Phi^t Phi dm  =  Sum Upsilon_kl(i,i)
Mgg=zeros(nf,nf);
if ~isempty(IW_U)
    for i=1:nf; for j=1:nf
        Mgg(i,j)=sum( (IW_U{i}(1,:).*U{j}(1,:) + IW_U{i}(2,:).*U{j}(2,:) + IW_U{i}(3,:).*U{j}(3,:)  ).*m   );
    end; end;
else
    for i=1:nf; for j=1:nf
        Mgg(i,j)=trapzs(s_span, (U{i}(1,:).*U{j}(1,:) + U{i}(2,:).*U{j}(2,:) + U{i}(3,:).*U{j}(3,:) ).*m);
    end; end;
end
% Adding torsion contribution if any
Mgg=Mgg+diag(GMJxx); 
% Ensuring orthogonality if requested
if bOrth
    Mgg=Mgg.*eye(nf);
end






% --- Building M and making it symmetric
MM=zeros(6+nf, 6+nf);
MM(1:3      , 1:3)      = Mass*eye(3); % Mxx
MM(1:3      , 4:6)      = Mxt;
MM(1:3      , 6+(1:nf)) = Mxg;
MM(4:6      , 4:6)      = Mtt;
MM(4:6      , 6+(1:nf)) = Mtg;
MM(6+(1:nf) , 6+(1:nf)) = Mgg;
MM=triu(MM)+triu(MM,1)';



function zz = trapzs(~,yy)
    % yy needs to be a line vector
    %zz =  (yy(1:nSpan-1) + yy(2:nSpan))/2 * dxx;
    zz = sum(yy.*IW);
end % trapz1d
% function zz = trapz1d(xx,yy)
%     xx=xx(:)';
%     yy=yy(:);
%     n = length(yy);
%     % Trapezoid sum computed with vector-matrix multiply.
%     % NOTE: x and y need have different dimension below
%     zz = diff(xx) * (yy(1:n-1) + yy(2:n))/2;
% end % trapz1d

end
