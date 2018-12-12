function [MM] = fGMBeamShapeIntegrals(gzf,Psi,Upsilon_kl,Sigma_kl,sigma_kl,sigma,Mass,ImomX,I_Jxx,GM_Jxx,bOrth,bAxialCorr,M1)

% NOTATIONS:
%   r or x     : 3 first translation dof
%   t or theta : 3 Rotational DOFS
%   f or g     : 3 flexible dofs
% 
% AUTHOR: E. Branlard

nf=size(Psi,2); % number of modes

% --- Optional argument
if ~exist('ImomX','var');  ImomX=0; end;
if ~exist('I_Jxx','var');  I_Jxx=zeros(1,nf); end;
if ~exist('GM_Jxx','var'); GM_Jxx=zeros(1,nf); end;
if ~exist('bOrth','var');  bOrth=false; end;
if ~exist('bAxialCorr','var');  bAxialCorr=false; end;



%%
gzf=gzf(:);


%% --- Mass Matrix 

% --- Mxx
Mxx= Mass*eye(3);

% --- Mxg = S    -  M x-g =  Psi
Mxg=Psi;
if bAxialCorr
    if ~exist('M1','var');  error('Provide M1 for axial corr'); end;
    for j=1:nf
        % KEEP ME
        %m16j  = GMVa7p(j)+GMVa7(j,:)*gzf;
        %m16j(j) = trapz(s_span, PhiU{j}(2,:).*peqy_tot  + PhiU{j}(3,:).*peqz_tot);
        %m16j(j) = trapz(s_span, -V{j}(2,:).*V_tot(2,:).*FXG  - V{j}(3,:).*V_tot(3,:).*FXG); 
        Mxg(1,j)= M1.m17p(j)+M1.m17g(j,:)*gzf;
    end
end

% --- Mxt =-Skew(I1+Sg)      -     M x-theta =-Skew(sigma+Psi g)
% Sxt=  -I1+S(:,1:nf)*gzf(1:nf);
% Mxt =- [0 -Sxt(3) Sxt(2) ; Sxt(3) 0 -Sxt(1) ; -Sxt(2) Sxt(1) 0 ]; % 
Sxt=sigma+Psi(:,1:nf)*gzf(1:nf);
Mxt = - [0 -Sxt(3) Sxt(2) ; Sxt(3) 0 -Sxt(1) ; -Sxt(2) Sxt(1) 0 ]; % - Skew matrix

if bAxialCorr
    % KEEP ME
    % m15 = -GMVa5p - GMVa5*gzf;
    % m16 =  GMVa6p + GMVa6*gzf;
    % m15 = trapz(s_span,  -s_span .* peqz_tot);
    % m16 = trapz( s_span,  s_span .* peqy_tot);
    %m15=-trapz(-s_span, V_tot(3,:).*FXG);
    %m16=-trapz( s_span, V_tot(2,:).*FXG);
    Mxt(1,2)= -M1.m15p - M1.m15g*gzf;
    Mxt(1,3)=  M1.m16p + M1.m16g*gzf;
end





% --- Mgg
Mgg= Upsilon_kl{1,1}+Upsilon_kl{2,2}+Upsilon_kl{3,3};
% Adding torsion contribution
Mgg=Mgg+diag(GM_Jxx); 
if bOrth % Keeping only diagonal elements
    Mgg=Mgg.*eye(nf);
end
% Adding Mtop : TODO not general!
% Mgg=Mgg+eye(nf)*Mtop;
% --- M t-t
% KEEP ME - Mtt  Brute force
%for i=1:length(o.s_span);
%    s=o.s_P0(:,i);
%    for j=1:o.nf
%        s=s+ o.gzf(j) * o.PhiU{j}(:,i);
%    end
%    DM(i,1:3,1:3)=-fSkew( s ) * fSkew(s) *o.m(i);
%end
%for i=1:3
%    for j=1:3
%        Mtt(i,j)=trapz(o.s_span, squeeze(DM(:,i,j)));
%    end
%end
Mtt=zeros(3,3);
Mx= sigma_kl(1,1) + 2*Sigma_kl{1,1}*gzf + gzf'*Upsilon_kl{1,1}*gzf ;
My= sigma_kl(2,2) + 2*Sigma_kl{2,2}*gzf + gzf'*Upsilon_kl{2,2}*gzf ;
Mz= sigma_kl(3,3) + 2*Sigma_kl{3,3}*gzf + gzf'*Upsilon_kl{3,3}*gzf ;
Mtt(1,1) = My + Mz + ImomX;
Mtt(2,2) = Mx + Mz ;
Mtt(3,3) = Mx + My ;
Mtt(1,2) = -(sigma_kl(2,1) + (Sigma_kl{1,2}+Sigma_kl{2,1})*gzf + gzf'*Upsilon_kl{2,1}*gzf);
Mtt(1,3) = -(sigma_kl(3,1) + (Sigma_kl{1,3}+Sigma_kl{3,1})*gzf + gzf'*Upsilon_kl{3,1}*gzf);
Mtt(2,3) = -(sigma_kl(2,3) + (Sigma_kl{2,3}+Sigma_kl{3,2})*gzf + gzf'*Upsilon_kl{3,2}*gzf);
Mtt=triu(Mtt)+triu(Mtt,1)';
% --- Mtg
Mtg=zeros(3,nf);
Mtg(1,1:nf)= gzf'*(Upsilon_kl{2,3}-transpose(Upsilon_kl{2,3})) + (Sigma_kl{2,3}-Sigma_kl{3,2}); % NOTE: primed required for Sigma_kl
Mtg(2,1:nf)= gzf'*(Upsilon_kl{3,1}-transpose(Upsilon_kl{3,1})) + (Sigma_kl{3,1}-Sigma_kl{1,3});
Mtg(3,1:nf)= gzf'*(Upsilon_kl{1,2}-transpose(Upsilon_kl{1,2})) + (Sigma_kl{1,2}-Sigma_kl{2,1});
% Adding torsion contribution if any
Mtg(1,1:nf)=Mtg(1,1:nf)+I_Jxx(1:nf);

% --- Building M and making it symmetric
MM=zeros(6+nf, 6+nf);
MM(1:3      , 1:3)      = Mxx;
MM(1:3      , 4:6)      = Mxt;
MM(1:3      , 6+(1:nf)) = Mxg;
MM(4:6      , 4:6)      = Mtt;
MM(4:6      , 6+(1:nf)) = Mtg;
MM(6+(1:nf) , 6+(1:nf)) = Mgg;
MM=triu(MM)+triu(MM,1)';

