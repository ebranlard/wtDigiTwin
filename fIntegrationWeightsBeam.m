function [IW,IW_x,IW_m,IW_xm,IW_U,Mass,Smom,Smom1,Imom] = fIntegrationWeightsBeam(s_span,m,PhiU,PhiV,PhiK)
% Returns integration weights convenient to integrate functions along the span of the beam
%
% 
% NOTE: Smom and Smom1 should be equal in theory
%
% AUTHOR: E.Branlard

% Ensuring line vectors
s_span = s_span(:)';
m      = m(:)';

% Usefull variables
nf    = length(PhiU)  ;
nSpan = length(s_span);

% Making sure modes have the proper dimensions
if size(PhiU{1},1)~=3     ;  error('The cells in PhiU need to be 3 x n');end;
if size(PhiU{1},2)~=nSpan ;  error('The dimension of PhiU needs to consistent with s_span');end;
if length(m)~=nSpan       ;  error('The length of m needs to be consistent with s_span');end;

% --- Span Integration weights  IW and IW_x 
% Assuming a fonction f that varies linearly
%     IW   is such that \int   f(x) dx = \Sum IW  [i] F[i]
%     IW_x is such that \int x.f(x) dx = \Sum IW_x[i] F[i]
% The equations are written such that s_span(1) is not necessary 0
IW  =zeros(1,nSpan); 
IW_x=zeros(1,nSpan); 
for I=1:nSpan-1 
    L           = s_span(I+1) - s_span(I);
    IW  (I)     = IW(I) + L/2          ;
    IW  (I+1)   = L/2                        ;
    IW_x(I)   = IW_x(I) + (s_span(I)/2 + L/6)*L;
    IW_x(I+1) = (s_span(I)/2 + L/3)*L;
end;
% --- Mass Integration weights IW_m IW_xm - CAN PROBABLY BE IMPROVED
%     IW_m  is such that \int   f(x).m(x) dx = \Sum IW_m [i] F[i]
%     IW_xm is such that \int x.f(x).m(x) dx = \Sum IW_xm[i] F[i]
IW_m=zeros(1,nSpan);
IW_xm=zeros(1,nSpan); 
for I=1:nSpan
    IW_m(I)  = IW(I)  *m(I);
    IW_xm(I) = IW_x(I)*m(I);
end;
% --- Shape function Integration Weight  
%    IW_U is such that  \int f(x).Ue1(x) dx = \Sum GFW[1,i] f[i]
IW_U=cell(1,nf);
for j=1:nf
    W=zeros(3,nSpan);
    for k=1:3
        U=PhiU{j}(k,:);
        V=PhiV{j}(k,:);
        K=PhiK{j}(k,:);
        for i=1:nSpan-1
          L = s_span(i+1)-s_span(i);
          if k==1
              W(1,i)   =  W(1,i)+(V(i) / 2 + (1 * K(i)/ 8   + K(i+1)/24)*L)*L;
              W(1,i+1) =         (V(i) / 2 + (5 * K(i)/ 24  + K(i+1)/8 )*L)*L;
          else
              W(k,i)   =  W(k,i)+(U(i)/2+(V(i)/6+(K(i)/30    +K(i+1)/120)*L)*L)*L;
              W(k,i+1) =         (U(i)/2+(V(i)/3+(K(i)*11/120+K(i+1)/30 )*L)*L)*L;
          end
        end; % loop on i
    end; %loop on k
    IW_U{j}=W;
end

% Calculation of mass, static moment and moment of inertia around X = 0 (RC)
Mass = 0;
Smom = 0;
Imom = 0;
Smom1 = 0;
for I=1:nSpan
    Mass  = Mass  + IW_m (I)          ; % Mass = \int     m(x) dx
    Smom  = Smom  + IW_xm(I)          ; % Smom = \int x  .m(x) dx
    Smom1 = Smom1 + IW_m (I)*s_span(I); % Smom1= \int x  .m(x) dx
    Imom  = Imom  + IW_xm(I)*s_span(I); % Imom = \int x^2.m(x) dx
end
