function out= problem6(x,mode)

n = length(x);
if mode==1               % compute F(x)

   A = gallery('tridiag', n, 0, 2.5, -1);
    Fx = A*x +0.5*sin(x) - 1;
  out=Fx;
elseif  mode==2          % compute the projection
    out=max(x,0);
end
end    
        
    
    


