function [Psi,Upsilon_kl,Sigma_kl,sigma_kl,sigma,Mass,Jxx,I_Jxx,GMJxx,M1] = fShapeIntegralsBeam(s_P0_in,s_span_in,m_in,jxxG_in,PhiU_in,PhiV_in,varargin)
    % NOTE: 
    %  - s_P0: location of the center of gravity, all the mass m0 is concentrated there
    %  -  m0  is mass per length, as function of the spanwise coordinate, assumed to be x
    %  -  jxxG is the add moment of inertia per length around the center of gravity
    %
    % AUTHOR: E.Branlard

    % Usefull variables
    nSpan = length(s_span_in);
    nf    = length(PhiU_in)  ;
    % --- Optional arguments
    %p=inputParser();
    p=fInputParser();
    p.KeepUnmatched=true;
    p.addParameter('IW_U',[]);
    p.addParameter('IW_xm',[]);
    p.addParameter('bInterp',false);
    p.addParameter('bAxialCorr',false);
    p.addParameter('V0',[]); % For Axial Corr
    p.parse(varargin{:}); p=p.Results;

    % --- emptying value if not provided
    if ~exist('PhiV_in','var'); PhiV_in=[]; end;
    % --- Default values if empty
    if isempty(s_P0_in); s_P0_in=zeros(3,nSpan); s_P0_in(1,:)=s_span_in; end;
    if isempty(jxxG_in); jxxG_in=zeros(1,nSpan); end;
    if isempty(PhiV_in) 
        PhiV_in=cell(1,nf);
        for j=1:nf
            PhiV_in{j}=zeros(3,nSpan);
        end
    end


    % Ensuring line vectors
    s_span_in = s_span_in(:)';
    m_in      = m_in(:)';

    % Making sure modes have the proper dimensions
    if nf>0
        if size(PhiU_in{1},1)~=3     ;  error('The cells in PhiU need to be 3 x n');end;
        if size(PhiU_in{1},2)~=nSpan ;  error('The dimension of PhiU needs to consistent with s_span');end;
    end
    if length(m_in)~=nSpan       ;  error('The length of m needs to be consistent with s_span');end;

    % --- Interpolation
    if p.bInterp
        nSpan=400;
        s_span = linspace(s_span_in(1),s_span_in(end),nSpan);
        m      = interp1(s_span_in,m_in,s_span)             ;
        jxxG   = interp1(s_span_in,jxxG_in,s_span)             ;

        s_P0      = zeros(3,nSpan);
        s_P0(1,:) = interp1(s_span_in,s_P0_in(1,:),s_span);
        s_P0(2,:) = interp1(s_span_in,s_P0_in(2,:),s_span);
        s_P0(3,:) = interp1(s_span_in,s_P0_in(3,:),s_span);

        PhiU=cell(1,nf);
        for j=1:nf
            PhiU{j}=zeros(3,nSpan);
            PhiU{j}(1,:) = interp1(s_span_in,PhiU_in{j}(1,:),s_span);
            PhiU{j}(2,:) = interp1(s_span_in,PhiU_in{j}(2,:),s_span);
            PhiU{j}(3,:) = interp1(s_span_in,PhiU_in{j}(3,:),s_span);
        end
        PhiV=cell(1,nf);
        for j=1:nf
            PhiV{j}=zeros(3,nSpan);
            PhiV{j}(1,:) = interp1(s_span_in,PhiV_in{j}(1,:),s_span);
            PhiV{j}(2,:) = interp1(s_span_in,PhiV_in{j}(2,:),s_span);
            PhiV{j}(3,:) = interp1(s_span_in,PhiV_in{j}(3,:),s_span);
        end
        if ~isempty(p.V0)
            V0_in=p.V0;
            p.V0(1,:) = interp1(s_span_in,V0_in(1,:),s_span);
            p.V0(2,:) = interp1(s_span_in,V0_in(2,:),s_span);
            p.V0(3,:) = interp1(s_span_in,V0_in(3,:),s_span);
        end
        if ~isempty(p.IW_U); error('Cannot use interp and IW_U'); end
        if ~isempty(p.IW_xm); error('Cannot use interp and IW_xm'); end
    else
        jxxG   = jxxG_in;
        s_span = s_span_in;
        m      = m_in     ;
        PhiU   = PhiU_in  ;
        PhiV   = PhiV_in  ;
        s_P0   = s_P0_in  ;
    end % interpolation





    % Mass integration
    Mass  = trapz(s_span, m);

    Jxx = trapz(s_span, jxxG);
    GMJxx=zeros(1,nf);
    I_Jxx=zeros(1,nf);
    for j=1:length(PhiV) % NOTE: only computed if V not empty!
        GMJxx(j)=trapz(s_span,PhiV{j}(1,:).*jxxG.*PhiV{j}(1,:));
        %I_Jxx(j)=sum(Bld.GFTJ.*Bld.VEVx(j,:)); 
        I_Jxx(j)=trapz(s_span,jxxG.*PhiV{j}(1,:));
    end
    % Converting dm into a matrix for convenience
    dm=repmat(m,3,1);
    % Shape Integral I1 = \int s_P0 dm - sigma
    sigma=trapz(s_span , s_P0 .* dm  ,  2 );
    % Shape Integral S - Psi
    Psi=NaN(3,nf);
    for j=1:nf
        Psi(:,j) = trapz(s_span, PhiU{j} .* dm , 2);
    end
    % --- KEEP ME "match my equations"
%         % Building the Matrix of ShapeFunctions
%         SMat = zeros(nSpan,3,nf);
%         for i=1:nSpan
%             for j=1:nf
%                 SMat(i,1:3,j)=PhiU{j}(:,i);
%             end
%         end
%         Ssij = zeros(nSpan,3,3,nf,nf);
%         for i=1:nSpan
%             for k=1:3
%                 Sk=squeeze(SMat(i,k,:))';
%                 for l=1:3
%                     Sl=squeeze(SMat(i,l,:))';
%                     Ssij(i,k,l,:,:)=Sk'*Sl; 
%                 end
%             end
%         end
    % --- Skl - Shape Integral -  Upsilon_kl
    Upsilon_kl=cell(3,3);
    for k=1:3
        for l=1:3
            Upsilon_kl{k,l}=zeros(nf,nf);
            for i=1:nf
                for j=1:nf
                    %Skl(k,l,i,j)=trapz(s_span, squeeze(Ssij(:,k,l,i,j))'.*m );      % NOTE KEEP ME
                    if ~isempty(p.IW_U)
                        Upsilon_kl{k,l}(i,j)=sum  (      p.IW_U{i}(k,:).*PhiU{j}(l,:).*m   );
                    else
                        Upsilon_kl{k,l}(i,j)=trapz(s_span, PhiU{i}(k,:).*PhiU{j}(l,:).*m   );
                    end
                end
            end
        end
%         if p.bOrth 
%             Upsilon_kl(k,k,:,:)=Upsilon_kl(k,k,:,:).*eye(3,3);
%         end
    end
    % --- Ikl - Shape Integral  - Sigma_kl
    Sigma_kl=cell(3,3);
    for k=1:3 % for beam only non zero for k=1
        for l=1:3
            Sigma_kl{k,l}=zeros(1,nf); % NOTE line vector
            for j=1:nf
                if ~isempty(p.IW_xm) && k==1
                    Sigma_kl{k,l}(j)=sum(            p.IW_xm .* PhiU{j}(l,:)     );
                else
                    Sigma_kl{k,l}(j)=trapz(s_span, s_P0(k,:) .* PhiU{j}(l,:) .* m);
                end
            end
        end
    end
    % --- ikl - Shape Integral  - sigma_kl
    sigma_kl=zeros(3,3);
    for k=1:3
        for l=1:3
            if ~isempty(p.IW_xm) && k==1
                sigma_kl(k,l)=sum(            p.IW_xm .* s_P0(l,:)     );
            else
                sigma_kl(k,l)=trapz(s_span, s_P0(k,:) .* s_P0(l,:) .* m);
            end
        end
    end


    % --- Axial correction (line 1 or mass matrix)
    if p.bAxialCorr
        if nf>0
            if isempty(p.V0);  error('Provide V0 for axial corr'); end;
        end
        % --- Variable required for line 1 (axial correction)
        % FXG from axial acceleration AX = 1 
        FXG=fcumtrapzlr(s_span,m);
        % %peqy_tot =  K_tot(2,:).*FXG - V_tot(2,:).*m;
        % %peqz_tot =  K_tot(3,:).*FXG - V_tot(3,:).*m;
        % % Equivalent force for unit axial acceleration and unit deflections
        % PACCy=zeros(nf,nSpan);
        % PACCz=zeros(nf,nSpan);
        % for j=1:nModes
        %   PACCy(j,:) = PhiK{j}(2,:).*FXG-PhiV{j}(2,:).*m;
        %   PACCz(j,:) = PhiK{j}(3,:).*FXG-PhiV{j}(3,:).*m;
        % end
        % % Prebend contribution
        % PACCyp=K0(2,:).*FXG - V0(2,:).*m;
        % PACCzp=K0(3,:).*FXG - V0(3,:).*m;

        M1.m15g = zeros(1,nf) ;
        M1.m16g = zeros(1,nf) ;
        M1.m17g = zeros(nf,nf);
        M1.m17p = zeros(1,nf) ;
        for j=1:nf
            for k=1:nf
                %GMVa7(J,k)=sum(Bld.GFWz(J,:).*Bld.PACCz(k,:)+Bld.GFWy(J,:).*Bld.PACCy(k,:));
                %GMVa7(j,k)=sum(IW_U{j}(3,:).*PACCz(k,:)+IW_U{j}(2,:).*PACCy(k,:));
                M1.m17g(j,k)=-trapz(s_span, (PhiV{j}(2,:).*PhiV{k}(2,:) + PhiV{j}(3,:).*PhiV{k}(3,:) ).*FXG);
            end
            %GMVa7p(j)= sum( Bld.GFWz(j,:).*Bld.PACCzp + Bld.GFWy(j,:).*Bld.PACCyp  ); 
            %GMVa7p(j)= sum(  IW_U{j}(3,:).*PACCzp + IW_U{j}(2,:).*PACCyp );
            M1.m17p(j)=-trapz(s_span, (PhiV{j}(2,:).*p.V0(2,:) + PhiV{j}(3,:).*p.V0(3,:) ).*FXG);
            %GMVa5(j) = sum(IW_x.*PACCz(j,:));
            M1.m15g(j)=-trapz(s_span, PhiV{j}(3,:).*FXG);
            M1.m16g(j)=-trapz(s_span, PhiV{j}(2,:).*FXG);
           % GMVa6(j) = sum(IW_x.*PACCy(j,:));
        end
        %GMVa5p=sum(IW_x.*PACCzp);
        %GMVa6p=sum(IW_x.*PACCyp);
        if nf>0
            M1.m15p=-trapz(s_span, p.V0(3,:).*FXG);
            M1.m16p=-trapz(s_span, p.V0(2,:).*FXG);
        else
            M1.m15p=0;
            M1.m16p=0;
        end
    else
        M1=[];
    end
end
