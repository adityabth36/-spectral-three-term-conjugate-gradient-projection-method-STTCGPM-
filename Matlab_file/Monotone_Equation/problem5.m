function out= problem5(x,mode)

 n = length(x);
if mode==1               % compute F(x)
  Fx=ones(n,1);
    Fx(1) = 2*x(1)+x(2) ;
  for i=2:n-1
      Fx(i) = cos(x(i-1))+2.5*x(i)-2+x(i+1) ;
  end
  Fx(n) = 2*x(n) -x(n-1);
  out=Fx;
elseif  mode==2          % compute the projection
    out=max(x,-1);
end  
end
    
        
    
    


