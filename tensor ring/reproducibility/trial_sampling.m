clc
clear
close all
%% input arguments
d=6;
I=6;
TRr=3;
sr=0.05:0.05:0.95;
tnsr=TR_rand(I*ones(1,d),TRr*ones(1,d));
num_trial=100;
mu_min=1e-12;
err_c=zeros(d/2,length(sr),num_trial);
%% counting
for l=1:d/2
    fprintf('l=%d\n',l);
    for i=1:length(sr)
        fprintf('i=%d\n',i);
        for j=1:num_trial
            if mod(j,20)==0
                fprintf('j=%d\n',j);
            end
            % sampling
            P=sampling_uniform(tnsr,sr(i));
            % solve problem via FPCA
            con=TRcompletion_FPCA_bound(tnsr,P,1,l,mu_min,false);
            err_c(l,i,j)=con;
        end
    end
end
%% visualize the results
err_s=sum(err_c,3)/num_trial;
figure(1);
plot(sr*100,err_s','linewidth',1.5);
xlabel('SR(%)');
ylabel('P');
%title('3*3*3*3*3*3*3*3');
legend('l=1','l=2','l=3','location','southeast');