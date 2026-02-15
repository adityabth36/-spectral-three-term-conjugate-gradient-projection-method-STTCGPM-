function out= problem4(x,mode)


n = length(x);
if mode==1
  Fx=ones(n,1);
  Fx(1)=2*x(1)+sin(x(1))-1;
  for i=2:n-1
    Fx(i)=2*x(i-1)+4*x(i)+sin(x(i))-1;
  end
  Fx(n)=2*x(n)+sin(x(n))-1;
  out=Fx;
elseif  mode==2
    out=max(x,0);
end
end
    
        
    
    


