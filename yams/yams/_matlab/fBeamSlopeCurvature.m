function [V,K]=fBeamSlopeCurvature(s_span,U,Vref,Kref,tol) 
    % INPUTS
    %   s_span : 1 x nSpan vector of spanwise position along the beam
    %   U      : cell(1 x nf) of (3 x nSpan) Beam deflections
    % 
    % OPTIONAL INPUTS
    %   If provided, these values are returned but checked against the numerical values
    %   Vref, Kref: like V and K. 
    %   tol: tolerance
    % OUPUTS:
    %   V : Slope, same dimensions as U
    %   K : Curvature, same dimensions as U
    %
    % AUTHOR: E. Branlard

    dx=s_span(2)-s_span(1);

    % Gradient implementation is for regular grid
    if length(unique(round(diff(s_span)*1000)))>1; 
        if bAbortAllowed
            error('Gradient below is for a regular grid'); 
        end
    end

    bDontReturnCell=false;
    if ~iscell(U)
        bDontReturnCell=true;
        U={U}; % temporarly comverting to cell 
        if nargin>=3 && ~isempty(Vref)
            Vref={Vref}; % temporarly comverting to cell 
            Kref={Kref}; % temporarly comverting to cell 
        end
    end

    nf=length(U);
    V = cell(1,nf); 
    K = cell(1,nf);
    % Computing V and K
    for j=1:nf
        V{j} = zeros(size(U{j}));
        K{j} = zeros(size(U{j}));
        for k=2:3 % start at 2, important we cant compute the torsional/twist angle with the mean line
            V{j}(k,:) = fgradient_regular(U{j}(k,:),4,dx);
            K{j}(k,:) = fgradient_regular(V{j}(k,:),4,dx);  % TODO improve me
        end
    end
    % 
    if nargin>=3  && ~isempty(Vref)
        % VERY IMPORTANT!!!!!!!!!!!!!!!!! COPY V0
        for j=1:nf
            V{j}(1,:)=Vref{j}(1,:);
            K{j}(1,:)=Kref{j}(1,:);
        end
        if isempty(Vref)
            fprintf('[WARN] Using numerically computed slope and curvature (V and K)\n');
        else % performing a sanity check
            for j=1:nf
                if norm(Vref{j}-V{j}) > tol; fprintf('[WARN] Slope V of mode %d seems to be off: %f\n',j,err); end
                if norm(Kref{j}-K{j}) > tol; fprintf('[WARN] Curv  K of mode %d seems to be off: %f\n',j,err); end
            end
            % returning the inputs...
            V=Vref;
            K=Kref;
        end
    end

    if bDontReturnCell
        V=V{1};
        K=K{1};
    end

