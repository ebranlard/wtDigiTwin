function B=fCreateBodyUniformBeam(Name,nShapes,nSpan,L,EI,m,Mtop,varargin)
% 
% 
% 
p = fInputParser();
% --- Key, value parameters
p.addParameter('bCompatibility',false);
p.addParameter('bAxialCorr',true);
p.addParameter('bStiffnessFromGM',false);
p.addParameter('bUseShapeIntegrals',true);
p.addParameter('jxx',[],@isnumeric);
p.addParameter('GKt',[],@isnumeric);
p.parse(varargin{:});
p=p.Results;



x=linspace(0,L,nSpan);
A=1;
rho=A*m;

% --------------------------------------------------------------------------------}
%% --- Determination of modes (theory)
% --------------------------------------------------------------------------------{
% --- Using external function
% [freq,~,U,V,K] = fUniformBeamTheory('transverse-unloaded-clamped-free',EI,m,1,L,'x',x,'norm','tip_norm');

% --- Copy Paste external function
x0=x/L;
% case 'transverse-unloaded-topmass-clamped-free'
% The geometrical stiffning is not accounted for here
M    = rho *A * L;
B=zeros(1,nShapes); % Beta
for i=1:nShapes
    B(i)=fzero(@(x) 1+cosh(x)*cos(x) - x*Mtop/M *(sin(x)*cosh(x)-cos(x)*sinh(x)),(2*i-1)*pi/2, optimset('TolFun',1e-21));
end
SS  = sinh(B)+sin(B);
CC  = cosh(B)+cos(B); 
% Frequency
freq = (B/L).^2/(2*pi)*sqrt(EI /(rho*A)); 
% Mode shapes - Normalized and in dimension space
U=zeros(nShapes,length(x0)); V=zeros(nShapes,length(x0)); K=zeros(nShapes,length(x0));
for i=1:nShapes
    U(i,:) =             SS(i)*(cosh(B(i)*x0)-cos(B(i)*x0)) - CC(i)*(sinh(B(i)*x0)-sin(B(i)*x0)) ;
    fact=1/U(i,end);
    U(i,:) = U(i,:)*fact;
    V(i,:) = B(i)  *(SS(i)*(sinh(B(i)*x0)+sin(B(i)*x0)) - CC(i)*(cosh(B(i)*x0)-cos(B(i)*x0)))/L  *fact;
    K(i,:) = B(i)^2*(SS(i)*(cosh(B(i)*x0)+cos(B(i)*x0)) - CC(i)*(sinh(B(i)*x0)+sin(B(i)*x0)))/L^2*fact;
end

% --------------------------------------------------------------------------------}
%% ---  Setting up YAMS Body
% --------------------------------------------------------------------------------{
B=Body('Name',Name,'Type','FlexibleBeam');
B.s_span=x;
B.ModeOmegas=zeros(1,nShapes);
for j=1:nShapes
    B.PhiU{j}=zeros(3,nSpan); B.PhiV{j}=zeros(3,nSpan); B.PhiK{j}=zeros(3,nSpan);
    B.PhiU{j}(3,:) = U(j,:);
    B.PhiV{j}(3,:) = V(j,:);
    B.PhiK{j}(3,:) = K(j,:);
    B.ModeOmegas(j) = freq(j) * 2*pi;
end
B.checkModesOrthogonality();
% ---
B.m   =  m*ones(size(B.s_span));
B.EIy = EI*ones(size(B.s_span));
if ~isempty(p.jxx)
    B.jxxG=p.jxx*ones(size(B.s_span));
end
if ~isempty(p.GKt)
    B.GKt=p.GKt*ones(size(B.s_span));
end
% B.s_C0_inB = [L; 0; 0];
% B.bOrthogonal=bOrthogonal;
% B.DampingMethod='StiffnessProportional';
% B.DampingParams=0.50;
B.bAxialCorr         = p.bAxialCorr        ;
B.bStiffnessFromGM   = p.bStiffnessFromGM  ;
B.setCompatible(p.bCompatibility);
B.bUseShapeIntegrals   = p.bUseShapeIntegrals  ;
p.bUseShapeIntegrals
B.Mtop=Mtop;
B.init();
