function [x,RC,RE,run_time]=TR_ADMM_bound(tnsr,P,mu,flag)
%% initialize parameters
maxiter=500;
epsilon_x=1e-8;
epsilon_err=1e-4;
N=ndims(tnsr);
J=size(tnsr);
F0=norm(tnsr(:),2);
x=P.*tnsr;
x0=x;
L=ceil(N/2);
m=cell(L,1);
y=cell(L,1);
sk=zeros(L,1);
for n=1:L
    m{n}=x;
    y{n}=zeros(J);
    order=[n:N 1:n-1];
    sk(n)=min([prod(J(order(1:L))),prod(J(order(L+1:N)))]);
end
m_cs=zeros(J);
y_cs=zeros(J);
idx_o=logical(P);
RC=nan(maxiter,1);
RE=nan(maxiter,1);
%% ADMM algorithm
[FR,r]=evaluate_fr(x,P);
if FR>=3
    epsilon=1e-2;
elseif FR>=2
    epsilon=1e-3;
else
    epsilon=1e-4;
end
w=harmmean(r)./r/length(r);
t=cputime;
% main loop
for i=1:maxiter
    % update m^(n)
    for n=1:L
        order=[n:N 1:n-1];
        Z_temp=permute(x-y{n}/mu,order);
        Z=reshape(Z_temp,prod(J(order(1:L))),[]);
        [M,sk(n)]=shrink_matrix(Z,w(n)/mu,sk(n),epsilon,true);
        M_temp=reshape(M,J(order));
        m{n}=ipermute(M_temp,order);
        m_cs=m_cs+m{n};
        y_cs=y_cs+y{n};
    end
    % update x
    x=(m_cs+y_cs/mu)/L;
    x(idx_o)=tnsr(idx_o);
    m_cs=zeros(J);
    y_cs=zeros(J);
    % update y^(n)
    for n=1:L
        y{n}=y{n}+mu*(m{n}-x);
    end
    % evaluate recovery accuracy
    RC(i)=norm(x(:)-x0(:),2)/norm(x0(:),2);
    RE(i)=norm(x(:)-tnsr(:),2)/F0;
    if flag && mod(i,10)==0
        fprintf('Iteration=%d\tRC=%f\tRE=%f\n',i,RC(i),RE(i));
    end
    if RC(i)<epsilon_x || RE(i)<epsilon_err || RE(i)>10
        break
    end
    x0=x;
    mu=min(mu*1.028,1e10);% 1.028 for images, 1.05 for videos
end
run_time=cputime-t;
fprintf('running time=%fs\n',run_time);
end