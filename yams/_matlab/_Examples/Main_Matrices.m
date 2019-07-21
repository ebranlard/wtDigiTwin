%%
% NOTE: This script was used to setup the unittests for flexible beam in Python
clear all; close all; clc;
format short;
addpath('../');


nShapes_bld=3;
nSpan_bld=30;
L_bld = 60   ; EI_bld = 2E+10; m_bld = 5E+2 ;
GKt_bld=7e11; % [Nm2]
jxx_bld=1e5;  % [kg.m]


% --- Test 1 - Straight beam
Bld=fCreateBodyUniformBeam('Bld',nShapes_bld,nSpan_bld,L_bld,EI_bld,m_bld,0,'jxx',jxx_bld,'GKt',GKt_bld);
Bld.bUseShapeIntegrals=false;
Bld.setCompatible(false);
Bld.computeMassMatrix();
Bld.MM


% --- Test 2 - curved beam with V_tot
Bld=fCreateBodyUniformBeam('Bld',nShapes_bld,nSpan_bld,L_bld,EI_bld,m_bld,0,'jxx',jxx_bld,'GKt',GKt_bld);
Bld.bUseShapeIntegrals=false;
Bld.setCompatible(false);
Bld.s_G(2,:)= Bld.s_G(1,:)/20;
Bld.s_G(3,:)= Bld.s_G(1,:)/10;
Bld.V_tot=Bld.PhiV{1}
Bld.bAxialCorr=true;
Bld.computeMassMatrix();


% --- Test 3
% Bld=fCreateBodyUniformBeam('Bld',nShapes_bld,nSpan_bld,L_bld,EI_bld,m_bld,0,'jxx',jxx_bld,'GKt',GKt_bld);
% Bld.setCompatible(true);
% Bld.s_P0(2,:)= Bld.s_P0(1,:)/20;
% Bld.s_P0(3,:)= Bld.s_P0(1,:)/10;
% Bld.s_G0=Bld.s_P0;
% Bld.s_G=Bld.s_P0;
% disp('NO AXIAL CORR')
% Bld.bAxialCorr=false;
% Bld.computeMassMatrix();
% Bld.MM
% disp('WITH AXIAL CORR')
% Bld.bAxialCorr=true;
% Bld.computeMassMatrix();
% Bld.MM
% 
% disp('WITH V0 CORR')
% Bld.V_tot=Bld.PhiV{1}
% Bld.computeMassMatrix();
% Bld.MM
% Bld.KK
% 

% fGMBeamStandalone(s_G ,s_span,m,jxxG,PhiUG,PhiV)

% KK0 = fGKBeam(s_span,GKt,EIy,EIz,PhiK,PhiV,V0,K0);
