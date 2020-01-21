clc
clear
close all
%% input arguments
d=4;
I=20;
sr=0.7;
TRr=10%floor(0.2*I);
tnsr=TR_rand(I*ones(1,d),TRr*ones(1,d));
%% sampling
P=sampling_uniform(tnsr,sr);
%% solve problem via TR-ADMM/SiLRTC-TT
mu_min=TRr^2/norm(P(:).*tnsr(:),2);
[x,~,RE,run_time]=TR_ADMM(tnsr,P,mu_min,true);