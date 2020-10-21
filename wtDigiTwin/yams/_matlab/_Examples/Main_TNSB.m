%% 
clear all; close all; clc;
format shortE;
addpath('../');

%% Main Parameters
nSpan_twr  = 101;
nSpan_bld  = 61;
nShapes_twr = 1 ; % Max 12
nShapes_bld = 1 ; % Max 12

bSft=1;
bHubMass=1;
bBld=1;
bInit=1;     %1  <<No Init to give mass matrix
nB=2;            %2 

bUseShapeIntegrals=false;

nDOF=1+nShapes_twr+nShapes_bld*nB;
if bInit
    q=zeros(nDOF,1);
    q(1:nDOF)=1; % initial conditions
else
    q=zeros(nDOF,1);
end

%% --- Strucural and geometrical Inputs
z_KGnac = 2e+00;
z_RGhub      = 0    ;
z_RS    = 0            ;
z_SR    = -z_RS                ;
z_NS    = - 10;
z_NGnac = 2e+00;         ;
z_SGhub = z_SR+z_RGhub         ;

r_NGnac_inN = [0;0;z_NGnac]
r_NS_inN    = [0;0;z_NS]
r_SGhub_inS = [0;0;z_SGhub] 
r_SR_inS    = [0;0;z_SR]
r_RGhub_inS = -r_SR_inS+r_SGhub_inS

M_hub=1e5*bHubMass;
IR_hub=zeros(3,3); % [kg*m^2].
IR_hub(1,1) = 2e+05;
IR_hub(2,2) = 2e+05;
IR_hub(3,3) = 3e+05; 
IR_hub=IR_hub*bHubMass;

M_nac   = 4e+05;
I0_nac=zeros(3,3);
I0_nac(1,1)=7e+06;
I0_nac(2,2)=3e+06;
I0_nac(3,3)=1e+06;

L_twr = 100  ; EI_twr = 2E+12; m_twr = 9E+3 ;
L_bld = 60   ; EI_bld = 2E+10; m_bld = 5E+2 ;
GKt_bld=7e11; % [Nm2]
jxx_bld=1e5;  % [kg.m]

% Inertias not at COG...
IG_hub = fTranslateInertiaMatrix(IR_hub, M_hub, [0;0;0], r_RGhub_inS);
IG_nac=fTranslateInertiaMatrixToCOG(I0_nac, M_nac,-r_NGnac_inN );

%% Derived parameters 
iPsi=nShapes_twr+1;


% --------------------------------------------------------------------------------}
%% --- Creating bodies
% --------------------------------------------------------------------------------{
%% Blade Body
bCompat =false;
Bld1=fCreateBodyUniformBeam('Bld1',nShapes_bld,nSpan_bld,L_bld,EI_bld,m_bld,0,'jxx',jxx_bld,'GKt',GKt_bld,'bCompatibility',bCompat,'bUseShapeIntegrals',bUseShapeIntegrals);
Bld2=fCreateBodyUniformBeam('Bld2',nShapes_bld,nSpan_bld,L_bld,EI_bld,m_bld,0,'jxx',jxx_bld,'GKt',GKt_bld,'bCompatibility',bCompat,'bUseShapeIntegrals',bUseShapeIntegrals);
Bld3=fCreateBodyUniformBeam('Bld3',nShapes_bld,nSpan_bld,L_bld,EI_bld,m_bld,0,'jxx',jxx_bld,'GKt',GKt_bld,'bCompatibility',bCompat,'bUseShapeIntegrals',bUseShapeIntegrals);

%% ShaftHub Body 
Sft=fCreateBodyRigid('ShaftHub',M_hub,IG_hub,r_SGhub_inS);

%% Nacelle Body
Nac=fCreateBodyRigid('Nacelle',M_nac,IG_nac,r_NGnac_inN);

%% Tower Body
if nB==2
    Mtop=(Bld1.Mass+Bld2.Mass)*bBld + Sft.Mass + Nac.Mass;
else
    Mtop=(Bld1.Mass+Bld2.Mass+Bld3.Mass)*bBld + Sft.Mass + Nac.Mass;
end
Twr=fCreateBodyUniformBeam('Tower',nShapes_twr,nSpan_twr,L_twr,EI_twr,m_twr,Mtop,'bCompatibility',false,'bAxialCorr',false,'bUseShapeIntegrals',bUseShapeIntegrals)
%% Ground
Grd=fCreateBodyGround('Ground');

% --------------------------------------------------------------------------------}
%% --- Automated connection of bodies and system matrix computation
% --------------------------------------------------------------------------------{
Grd.connectTo(Twr,'Point',[0;0;0],'Type','Rigid');
Twr.connectTo(Nac,'Point','LastPoint','Type','Rigid');
Nac.connectTo(Sft,'Point',r_NS_inN,'Type','SphericalJoint','JointRotations',{'z'},'Orientation',fRotz(pi));
Sft.connectTo(Bld1,'Point',r_SR_inS,'Type','Rigid','Orientation',fRotz(0));
if nB==2
    Sft.connectTo(Bld2,'Point',r_SR_inS,'Type','Rigid','Orientation',fRotz(pi));
else
    Sft.connectTo(Bld2,'Point',r_SR_inS,'Type','Rigid','Orientation',fRotz(-2*pi/3));
    Sft.connectTo(Bld3,'Point',r_SR_inS,'Type','Rigid','Orientation',fRotz( 2*pi/3));
end

% Setting DOF index for all bodies and connections 
nq=Grd.setupDOFIndex(0);
if nq~=nDOF
    error('Incompatibility')
end
qdot=zeros(nq,1);
qddot=zeros(nq,1);


%% Kinematics
Grd.updateChildrenKinematicsNonRecursive(q,qdot,qddot);
Twr.updateChildrenKinematicsNonRecursive(q,qdot,qddot);
Nac.updateChildrenKinematicsNonRecursive(q,qdot,qddot);
Sft.updateChildrenKinematicsNonRecursive(q,qdot,qddot);

%% Mass and stiffness matris
M=zeros(nq,nq);
M=Grd.getFullM(M);
K=zeros(nq,nq);
K=Grd.getFullK(K);

%% Eigenvalue analysis
[Q,Lambda]=eig(K,M);
Omega2=diag(Lambda);
[Omega2,Isort]=sort(Omega2);
Q=Q(:,Isort);
f_eva= sqrt(Omega2)/(2*pi);
for i=1:length(f_eva);
    fprintf('f%d = %.3f \n',i,f_eva(i))
end

%% Display to screen 
disp('Tower M')
disp(Twr.MM)
disp('Blade M')
disp(Bld1.MM)

M(abs(M)<1e-9)=0;
disp('M')
disp(M)
disp('K')
disp(K)









% --------------------------------------------------------------------------------}
%% --- "Manual connection"
% --------------------------------------------------------------------------------{
nf_T=Twr.nf;

% link E-T
R_ET=eye(3);
B_T=[];
B_T_inT=B_T; % B_T_inT=R_ET'*B_T;
BB_T_inT=[B_T_inT zeros(6,nf_T) ; zeros(nf_T, size(B_T_inT,2)) eye(nf_T)];

% Update of mass matrix
MM_T= BB_T_inT'*Twr.MM*BB_T_inT;
KK_T= BB_T_inT'*Twr.KK*BB_T_inT;
% CyT=zeros(1,Twr.nf);
if nShapes_twr==2
    CyT= - [Twr.PhiV{1}(3,end) Twr.PhiV{2}(3,end)];
else
    CyT= - Twr.PhiV{1}(3,end);
end

% ---------------------------------------------
% Link T-N
r_TN=[L_twr; 0; 0];

if nShapes_twr ==2 
    Bx_TN=[0 0; 0 0; 1 1]; % Bhat_x^N
    Bt_TN=[0 0; CyT(1) CyT(2); 0 0]; % Bhat_x^N
    alpha_y = CyT*q(1:2);
else
    Bx_TN=[0; 0; 1]; % Bhat_x^N
    Bt_TN=[0; CyT(1);  0]; % Bhat_x^N
    alpha_y = CyT*q(1);
end

R_TN=fRoty(alpha_y);

R_EN=R_ET*R_TN; % <<<<<< This is where I could insert some fixed tilt

B_N=fBMatRecursion(B_T,[Bx_TN; Bt_TN],R_ET,r_TN);
B_N_inN = [R_EN' * B_N(1:3,:);  R_EN' * B_N(4:6,:)];

nf_N=0;
BB_N_inN =[B_N_inN zeros(6,nf_N) ; zeros(nf_N, size(B_N_inN,2)) eye(nf_N)];
% Update of mass matrix
MM_N= BB_N_inN'*Nac.MM*BB_N_inN;
KK_N= BB_N_inN'*Nac.KK*BB_N_inN;

% ---------------------------------------------
% Link N-S
q_psi=q(iPsi); % <<<<<<<<<<<<<<<<<<<<<<< Needs update
R_NS = fRotz(q_psi+pi); % Adding pi here , blade down
R_ES=R_EN * R_NS;

r_NS=R_EN*r_NS_inN;

Bx_NS=[0;0;0];
Bt_NS=[0;0;1];

B_S=fBMatRecursion(B_N,[Bx_NS; Bt_NS],R_EN,r_NS);
B_S_inS = [R_ES' * B_S(1:3,:);  R_ES' * B_S(4:6,:)];

nf_S=Sft.nf;
BB_S_inS =[B_S_inS zeros(6,nf_S) ; zeros(nf_S, size(B_S_inS,2)) eye(nf_S)];

MM_S= BB_S_inS'*Sft.MM*BB_S_inS;
KK_S= BB_S_inS'*Sft.KK*BB_S_inS;

% ---------------------------------------------
% Link S-B1
if nB==2
    R_SB1 = fRotz(0*pi+0);
    R_SB2 = fRotz(0*pi+pi);
    R_SB3 = zeros(3,3);
else
    R_SB1 = fRotz(0*pi+0);
    R_SB2 = fRotz(0*pi-2*pi/3);
    R_SB3 = fRotz(0*pi+2*pi/3);
end
R_EB1 = R_ES*R_SB1;
R_EB2 = R_ES*R_SB2;
R_EB3 = R_ES*R_SB3;

r_SR=R_ES*r_SR_inS;

B_R     = fBMatRecursion(B_S,[],R_ES,r_SR);
B_R_bis = fBMatTranslate(B_S,r_SR);

B_B1_inB1 = [R_EB1' * B_R(1:3,:);  R_EB1' * B_R(4:6,:)];
B_B2_inB2 = [R_EB2' * B_R(1:3,:);  R_EB2' * B_R(4:6,:)];
B_B3_inB3 = [R_EB3' * B_R(1:3,:);  R_EB3' * B_R(4:6,:)];

nf_B1=Bld1.nf;
nf_B2=Bld2.nf;
if nB==2
    nf_B3=0;
else
    nf_B3=Bld3.nf;
end
I_B1 = [eye(nf_B1)   zeros(nf_B2) zeros(nf_B3)];
I_B2 = [zeros(nf_B1) eye(nf_B2)   zeros(nf_B3)];
I_B3 = [zeros(nf_B1) zeros(nf_B2)   eye(nf_B3)];
nf_B=nf_B1+nf_B2+nf_B3;
BB_B1_inB1 =[B_B1_inB1 zeros(6,nf_B) ; zeros(nf_B1, size(B_B1_inB1,2)) I_B1];
BB_B2_inB2 =[B_B2_inB2 zeros(6,nf_B) ; zeros(nf_B2, size(B_B2_inB2,2)) I_B2];
if nB==3
    BB_B3_inB3 =[B_B3_inB3 zeros(6,nf_B) ; zeros(nf_B3, size(B_B3_inB3,2)) I_B3];
end

MM_B1 = BB_B1_inB1'*Bld1.MM*BB_B1_inB1;
KK_B1 = BB_B1_inB1'*Bld1.KK*BB_B1_inB1;
MM_B2 = BB_B2_inB2'*Bld2.MM*BB_B2_inB2;
KK_B2 = BB_B2_inB2'*Bld2.KK*BB_B2_inB2;
if nB==3
    MM_B3 = BB_B3_inB3'*Bld3.MM*BB_B3_inB3;
    KK_B3 = BB_B3_inB3'*Bld3.KK*BB_B3_inB3;
else
    MM_B3 = MM_B2*0;
    KK_B3 = KK_B2*0;
end

% Final assembly
bBldMass=1;
MM=(MM_B1+MM_B2+MM_B3)*bBldMass;
MM(1:iPsi,1:iPsi)=MM(1:iPsi,1:iPsi)+MM_S;
MM(1:nShapes_twr,1:nShapes_twr)=MM(1:nShapes_twr,1:nShapes_twr)+MM_T+MM_N;


KK=(KK_B1+KK_B2+KK_B3);
KK(1:iPsi,1:iPsi)=KK(1:iPsi,1:iPsi)+KK_S;
KK(1:nShapes_twr,1:nShapes_twr)=KK(1:nShapes_twr,1:nShapes_twr)+KK_T+KK_N;

%% Display to screen
MM(abs(MM)<1e-9)=0;
disp('M ("manually" built)')
disp(MM)
disp('K ("manually" build)')
disp(KK)
Twr.MM
