function out= problem1(x,mode)

 n = length(x);
if mode==1               % compute F(x)
    % a=10;
    h=1/(n+1);
  Fx=ones(n,1);
    Fx(1) = 2*x(1)+0.5*h^2*(x(1)+h)^3-x(2) ;
  for i=2:n-1
      Fx(i) = 2*x(i)+0.5*h^2*(x(i)+h*i)^3-x(i-1)+x(i+1) ;
  end
  Fx(n) =2*x(n)+0.5*h^2*(x(n)+h*n)^3-x(n-1) ;
  out=Fx;
elseif  mode==2          % compute the projection
    out=max(x,0);
end 