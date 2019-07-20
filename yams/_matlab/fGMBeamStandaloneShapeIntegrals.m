function [MM,Psi,Upsilon_kl,Sigma_kl,sigma_kl,sigma,Mass,ImomX,I_Jxx,GMJxx,M1] = fGMBeamStandaloneShapeIntegrals(s_G0,s_span,m,JxxG,PhiU,PhiV,gzf,varargin)
% NOTE: Beam assumed to be along x for now (only because of JxxG)
% 
% INPUTS
%  - s_G   : 3 x nSpan , location of UNDEFLECTED COG
%  - s_span : span integration variable (e.g. s_G(1,:))
%
%
% OPTIONAL INPUTS:
%  - JxxG0, if omitted, assumed to be 0
%  - PhiU , if omitted, then rigid body (6x6) mass matrix is returned
%  - gzf  , if omitted, then Shape integrals are not
%
% --- Optional arguments
if ~exist('JxxG','var');  JxxG=0; end;
if ~exist('PhiU','var');  PhiU=[]; end
if ~exist('PhiV','var');  PhiV=[]; end
if ~exist('gzf','var');   gzf=zeros(length(PhiU),1); end;

% --- Optional arguments
p=fInputParser();
p.KeepUnmatched=true;
p.addParameter('bOrth',false);
p.addParameter('bAxialCorr',false);
p.addParameter('V_tot',[]);
p.addParameter('V0',[]);
p.parse(varargin{:});
p=p.Results;

% --- Computing shape integrals
[Psi,Upsilon_kl,Sigma_kl,sigma_kl,sigma,Mass,ImomX,I_Jxx,GMJxx,M1] = fShapeIntegralsBeam(s_G0,s_span,m,JxxG,PhiU,PhiV,varargin{:});
% --- Computing Mass matrix from shape integrals
MM = fGMBeamShapeIntegrals(gzf,Psi,Upsilon_kl,Sigma_kl,sigma_kl,sigma,Mass,ImomX,I_Jxx,GMJxx,p.bOrth,p.bAxialCorr,M1);


