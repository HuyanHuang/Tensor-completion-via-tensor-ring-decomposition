function [A,B,tnsr]=tensor_ring(g,n,I)
% g is the cores of tensor ring, which is represented by a cell that contains N elements,
% each element is a tensor of size Rn x In x Rn+1
% A*B will yield a matrix of size In x (In+1...IN*I1...In-1) 
% tnsr will be exactly the same as the origin 
gA=permute(g{n},[2 1 3]);
A=reshape(gA,I(n),[]);
if nargout>1
    N=length(g);
    order=[n+1:N 1:n-1];
    gB=g{order(1)};
    for i=2:N-1
%         [Rn,In,~]=size(gB);
%         [~,In1,Rn2]=size(g{order(i)});
%         gtemp=zeros(Rn,In*In1,Rn2);
%         for j=1:Rn
%             for k=1:Rn2
%                 temp=squeeze(gB(j,:,:))*squeeze(g{order(i)}(:,:,k));
%                 gtemp(j,:,k)=temp(:);
%             end
%         end
%         gB=gtemp;
        [Rn,In,Rn1]=size(gB);
        [~,In1,Rn2]=size(g{order(i)});
        gB=reshape(reshape(gB,[],Rn1)*reshape(g{order(i)},Rn1,[]),[Rn,In*In1,Rn2]);
    end
    gB=permute(gB,[3 1 2]);
    B=reshape(gB,[],prod(I)/I(n));
end
if nargout>2
    tnsr_temp=reshape(A*B,I([n order]));
    tnsr=ipermute(tnsr_temp,[n order]);
end
end