clc
clear
close all
%% input arguments
d=4;
I=20;
TRr=2;
TTr=2;
sr=0.1:0.05:0.9;
tnsr_TR=TR_rand(I*ones(1,d),TRr*ones(1,d));
tnsr_TT=TG_rand(I*ones(1,d),TTr,[1,d]);
mu_min=1e-1;
num_trial=10;
err_TR=zeros(length(sr),num_trial);
err_TT=zeros(length(sr),num_trial);
%% sampling
for i=1:length(sr)
    fprintf('i=%d\n',i);
    for j=1:num_trial
        fprintf('j=%d\n',j);
        P=sampling_uniform(tnsr_TR,sr(i));
        [~,ranktube] = SVD_MPS_Rank_Estimation(P.*tnsr_TT,10^-2.5);
        alpha=ranktube;
        gamma=0.01*alpha;
        %% solve problem via TR-ADMM/SiLRTC-TT
        [r,fr]=evaluate_rfr(tnsr,P);
        x=TR_ADMM(tnsr_TR,P,mu_min,false);
        err_TR(i,j)=norm(x(:)-tnsr_TR(:),2)/norm(tnsr_TR(:),2);
        x=SiLRTC_TT(P.*tnsr_TT,find(P),alpha,gamma,500,1e-5);
        err_TT(i,j)=norm(x(:)-tnsr_TT(:),2)/norm(tnsr_TT(:),2);
        %[x,err,run_time]=TRcompletion_ALS_Pro(tnsr,P,TRr,true);
    end
end
%% visualize the results
figure(1);
boxplot(log10(err_TR)',sr*100);
hold on
boxplot(log10(err_TT)',sr*100);
plot(log10(mean(err_TR,2)),'linewidth',1.5);
plot(log10(mean(err_TT,2)),'linewidth',1.5);
ylim([-16.5,0.5]);
xlabel('SR(%)');
ylabel('lg(RE)');
legend('TRBU','SiLRTC-TT','Location','best');