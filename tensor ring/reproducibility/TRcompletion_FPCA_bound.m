function con=TRcompletion_FPCA_bound(tnsr,P,n,l,mu_min,flag)
%% compute parameters
N=ndims(tnsr);
J=size(tnsr);
%% initialize parameters
order=[n:N 1:n-1];
P_temp=permute(P,order);
p=reshape(P_temp,prod(J(order(1:l))),[]);
T_temp=permute(tnsr,order);
T=reshape(T_temp,prod(J(order(1:l))),[]);
t=1.999;
mu_max=norm(p.*T,2)*exp(-4);
sk=min(size(p.*T));
mu=mu_max;
epsilon=1e-3;
xtol=1e-4;
errtol=1e-6;
x=zeros(size(T));
x_old=-ones(size(T));
maxiter=500;
err=zeros(maxiter,1);
err(1)=inf;
F0=norm(tnsr(:),2);
num_cvg=0;
num_dvg=0;
i=2;
con=0;
tol=1e-3;
%% FPCA algorithm
while mu>mu_min
    while norm(x-x_old,'fro')/norm(x_old,'fro')>xtol*sqrt(mu/mu_max)
        x_old=x;
        y=x-t*p.*(x-T);
        %[u,s,v]=svds(y,sk);
        [u,s,v]=svd(y,'econ');
        u=u(:,1:sk);
        s=s(1:sk,1:sk);
        v=v(:,1:sk);
        s_vec=diag(s)-t*mu;
        ind=find(s_vec>0,1,'last');
        sk=find(s_vec>=epsilon*s_vec(1),1,'last');
        s_shrink=diag(s_vec(1:ind));
        x=u(:,1:ind)*s_shrink*v(:,1:ind)';
        x=(1-p).*x+p.*T;
        if norm(x,'fro')>norm(x_old,'fro')
            %sk=sk+1;
        end
        err(i)=norm(x-T,'fro')/F0;
        if flag
            fprintf('Iteratoin=%d\tRE=%f\n',i-1,err(i));
        end
        if err(i)<tol
            num_cvg=20;
            con=1;
            break
        end
        if err(i)>10
            num_dvg=20;
            con=0;
            break
        end
        if err(i)>err(i-1)
            num_dvg=num_dvg+1;
        else
            if 1-err(i)/err(i-1)<errtol
                num_cvg=num_cvg+1;
            end
        end
        i=i+1;
        if num_dvg>=20||num_cvg>=20
            break    
        end
    end
    if num_dvg>=20||num_cvg>=20
        break    
    end
    mu=mu*exp(-1);
end
end