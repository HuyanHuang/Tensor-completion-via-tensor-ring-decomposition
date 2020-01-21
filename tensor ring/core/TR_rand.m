function [tnsr,g]=TR_rand(I,R)
g=TRg_rand(I,R);
[~,~,tnsr]=tensor_ring(g,1,I);
end