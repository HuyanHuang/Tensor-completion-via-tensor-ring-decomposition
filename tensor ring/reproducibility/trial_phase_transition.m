clc
clear
close all
%% input arguments
d=4;
I=20;
TRr=2:19;
SR=0.02:0.02:0.98;
num_trial=10;
% err1=zeros(length(TRr),length(sr),num_trial);
load err_phase
% err2=zeros(length(TRr),length(sr),num_trial);
%% main loops
for i=5:8%length(TRr)
    i
    tnsr=TR_rand(I*ones(1,d),TRr(i)*ones(1,d));
    for j=10:25%length(sr)
        j
        for k=1:num_trial
            P=sampling_uniform(tnsr,SR(j));
%             X=M_ADMM(tnsr,P,[]);
%             x=reshape(X,size(tnsr));
            [~,~,RE]=TR_ADMM_bound(tnsr,P,TRr(i)^2/norm(P(:).*tnsr(:),2),false);
            err1(i,j,k)=min(RE);
%             [~,~,RE]=TR_ALS(tnsr,P,TRr(i),false);
%             err2(i,j,k)=min(RE);
        end
    end
end
%% pre-processing
save err_phase err1
err1_c=squeeze(mean(err1,3,'omitnan'));
% err2_c=squeeze(mean(err2,3,'omitnan'));
%% visualize results
figure;
imagesc([TRr(1) TRr(end)],[SR(1) SR(end)],log10(err1_c'));
xlabel('TR-rank');
ylabel('SR');
colormap(flipud(gray));
caxis([-4 0]);
colorbar('Ticks',[-4 0],'TickLabels',{'1','0'},'Direction','reverse');
set(gca,'YDir','normal');
%%
d_m=(TRr.^2.*(2*I^(d/2)-TRr.^2))./(I^d*SR');
figure;
imagesc([TRr(1) TRr(end)],[SR(1) SR(end)],d_m);
xlabel('TR-rank');
ylabel('SR');
colormap(flipud(gray));
caxis([0 1]);
colorbar('Ticks',[0 1],'TickLabels',{'0','>1'},'Direction','reverse');
set(gca,'YDir','normal');
title('df_M/m');
%%
d_t=(d*I*TRr.^2-d*TRr.^2+1)./(I^d*SR');
figure;
imagesc([TRr(1) TRr(end)],[SR(1) SR(end)],d_t);
xlabel('TR-rank');
ylabel('SR');
colormap(flipud(gray));
caxis([0 1]);
colorbar('Ticks',[0 1],'TickLabels',{'0','>1'},'Direction','reverse');
set(gca,'YDir','normal');
title('df_{TR}/m');