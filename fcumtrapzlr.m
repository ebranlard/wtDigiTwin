function I=fcumtrapzlr(x,y)
% Left-right cumtrapz, computes I(xi) = \int_xi^xn y(x) dx
% 
% AUTHOR: E. Branlard
%
x=x(:)';
y=y(:)';

I = - fliplr( cumtrapz( fliplr(x), fliplr(y)) ) ;
I = reshape(I,size(x));
