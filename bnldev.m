%function res=bnldev(N,p)
% Returns binomially distributed random deviates =
%    number of successes in N(i) trials with probability p(i) per
%    trial
% Input
%     N = vector of trial sizes
%     p = vector of probabilities, length(p)=length(N)
% Output
%     res = vector of deviates
% Implementation of algorithm in Numerical Recipes in Fortran
% Works efficiently for large N, as opposed to built-in binornd
function res=bnldev(N,pp)
    res=zeros(1,length(N));
    for i=1:length(N)
        res(i)=bnldev_scalar(N(i),pp(i));
    end
end

function res=bnldev_scalar(N,pp)
p=min(pp,1-pp);
mu=N*p;
if (N<25)
    res=sum(rand(1,N)<p);
elseif (mu < 1)
    g=exp(-mu);
    t=1;
    for j=0:N
        t=t*rand(1);
        if (t<g) 
            break;
        end
    end
    res=j;
else
    g=gammaln(N+1);
    sq=sqrt(2*mu*(1-p));
    t=0;
    while (rand(1) > t)
        em=-1;
        while (em<0 || em>N+1)
            y=tan(pi*rand(1));
            em=sq*y+mu;
        end
        em=floor(em);
        t=1.2*sq*(1+y^2)*exp(g-gammaln(em+1)-gammaln(N-em+1)+em*log(p)...
            +(N-em)*log(1-p));
    end
    res=em;
end
if (pp > p) 
    res=N-res;
end
end