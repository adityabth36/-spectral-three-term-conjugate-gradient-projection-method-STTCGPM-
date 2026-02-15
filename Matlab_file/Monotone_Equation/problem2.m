function out= problem2(x,mode)


n = length(x);
if mode==1
    out=x-sin(abs((x)-3));
elseif mode==2
    if sum(x)<=n && min(x)>=-1
        out=x;
    else
        out=quadprog(speye(n),-x,ones(1,n),n,[],[],-ones(n,1));
        % out=quadprog(eye(n),-x,[ones(1,n);zeros(n-1,n)],[n;zeros(n-1,1)],[],[],-1*ones(n,1),[]);
    end
end
end
