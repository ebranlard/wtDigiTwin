function I=fcumtrapzlr(x,y)
% Left-right cumtrapz, computes I(xi) = \int_xi^xn y(x) dx
% 
% AUTHOR: E. Branlard
%
x=x(:)';
y=y(:)';

I = - fliplr( cumtrapz( fliplr(x), fliplr(y)) ) ;
I = reshape(I,size(x));
%%
% I=zeros(size(y));
% for i=length(x)-1:-1:1
%     dx=x(i+1)-x(i);
%     I(i)=I(i+1)+ (y(i)+y(i+1))*dx/2;
% end

