clc
clear
close all
%% input arguments
sr=0.1:0.1:0.5;
num_trial=1;
flag=false;
%% image tests
L=length(sr);
% err_psnr=cell(1,8);
% run_time=cell(1,8);
load err_psnr1_VDT
load run_time1_VDT
warning off
for i=1:8
    i
    switch i
        case 1
            img=imread('kodim04.png');
            img=double(img);
            for j=1:L
                j
                for k=1:num_trial
                    % sampling
                    P=sampling_uniform(img,sr(j));
                    % solve problem via TR-ADMM
                    [run_time,err_psnr]=test_method1(img,P,run_time,err_psnr,i,j,k,1,flag);
                    % solve problem via TR-ALS
                    if k==1
                        [run_time,err_psnr]=test_method1(img,P,run_time,err_psnr,i,j,k,2,flag);
                    end
                    % solve problem via SiLRTC-TT
                    [run_time,err_psnr]=test_method1(img,P,run_time,err_psnr,i,j,k,3,flag);
                    % solve problem via LRTC-TNN
                    [run_time,err_psnr]=test_method1(img,P,run_time,err_psnr,i,j,k,4,flag);
                    % solve problem via FBCP
                    [run_time,err_psnr]=test_method1(img,P,run_time,err_psnr,i,j,k,5,flag);
                    % solve problem via HaLRTC
                    [run_time,err_psnr]=test_method1(img,P,run_time,err_psnr,i,j,k,6,flag);
                    % solve problem via STTC
                    [run_time,err_psnr]=test_method1(img,P,run_time,err_psnr,i,j,k,7,flag);
                    % solve problem via TRNNM
                    [run_time,err_psnr]=test_method1(img,P,run_time,err_psnr,i,j,k,8,flag);
                end
            end
        case 2
            img=imread('peppers.bmp');
            img=double(img);
            for j=1:L
                j
                for k=1:num_trial
                    % sampling
                    P=sampling_uniform(img,sr(j));
                    % solve problem via TR-ADMM
                    [run_time,err_psnr]=test_method1(img,P,run_time,err_psnr,i,j,k,1,flag);
                    % solve problem via TR-ALS
                    if k==1
                        [run_time,err_psnr]=test_method1(img,P,run_time,err_psnr,i,j,k,2,flag);
                    end
                    % solve problem via SiLRTC-TT
                    [run_time,err_psnr]=test_method1(img,P,run_time,err_psnr,i,j,k,3,flag);
                    % solve problem via LRTC-TNN
                    [run_time,err_psnr]=test_method1(img,P,run_time,err_psnr,i,j,k,4,flag);
                    % solve problem via FBCP
                    [run_time,err_psnr]=test_method1(img,P,run_time,err_psnr,i,j,k,5,flag);
                    % solve problem via HaLRTC
                    [run_time,err_psnr]=test_method1(img,P,run_time,err_psnr,i,j,k,6,flag);
                    % solve problem via STTC
                    [run_time,err_psnr]=test_method1(img,P,run_time,err_psnr,i,j,k,7,flag);
                    % solve problem via TRNNM
                    [run_time,err_psnr]=test_method1(img,P,run_time,err_psnr,i,j,k,8,flag);
                end
            end
        case 3
            img=imread('sailboat.bmp');
            img=double(img);
            for j=1:L
                j
                for k=1:num_trial
                    % sampling
                    P=sampling_uniform(img,sr(j));
                    % solve problem via TR-ADMM
                    [run_time,err_psnr]=test_method1(img,P,run_time,err_psnr,i,j,k,1,flag);
                    % solve problem via TR-ALS
                    if k==1
                        [run_time,err_psnr]=test_method1(img,P,run_time,err_psnr,i,j,k,2,flag);
                    end
                    % solve problem via SiLRTC-TT
                    [run_time,err_psnr]=test_method1(img,P,run_time,err_psnr,i,j,k,3,flag);
                    % solve problem via LRTC-TNN
                    [run_time,err_psnr]=test_method1(img,P,run_time,err_psnr,i,j,k,4,flag);
                    % solve problem via FBCP
                    [run_time,err_psnr]=test_method1(img,P,run_time,err_psnr,i,j,k,5,flag);
                    % solve problem via HaLRTC
                    [run_time,err_psnr]=test_method1(img,P,run_time,err_psnr,i,j,k,6,flag);
                    % solve problem via STTC
                    [run_time,err_psnr]=test_method1(img,P,run_time,err_psnr,i,j,k,7,flag);
                    % solve problem via TRNNM
                    [run_time,err_psnr]=test_method1(img,P,run_time,err_psnr,i,j,k,8,flag);
                end
            end
        case 4
            img=imread('lena.bmp');
            img=double(img);
            for j=1:L
                j
                for k=1:num_trial
                    % sampling
                    P=sampling_uniform(img,sr(j));
                    % solve problem via TR-ADMM
                    [run_time,err_psnr]=test_method1(img,P,run_time,err_psnr,i,j,k,1,flag);
                    % solve problem via TR-ALS
                    if k==1
                        [run_time,err_psnr]=test_method1(img,P,run_time,err_psnr,i,j,k,2,flag);
                    end
                    % solve problem via SiLRTC-TT
                    [run_time,err_psnr]=test_method1(img,P,run_time,err_psnr,i,j,k,3,flag);
                    % solve problem via LRTC-TNN
                    [run_time,err_psnr]=test_method1(img,P,run_time,err_psnr,i,j,k,4,flag);
                    % solve problem via FBCP
                    [run_time,err_psnr]=test_method1(img,P,run_time,err_psnr,i,j,k,5,flag);
                    % solve problem via HaLRTC
                    [run_time,err_psnr]=test_method1(img,P,run_time,err_psnr,i,j,k,6,flag);
                    % solve problem via STTC
                    [run_time,err_psnr]=test_method1(img,P,run_time,err_psnr,i,j,k,7,flag);
                    % solve problem via TRNNM
                    [run_time,err_psnr]=test_method1(img,P,run_time,err_psnr,i,j,k,8,flag);
                end
            end
        case 5
            img=imread('barbara.bmp');
            img=double(img);
            for j=1:L
                j
                for k=1:num_trial
                    % sampling
                    P=sampling_uniform(img,sr(j));
                    % solve problem via TR-ADMM
                    [run_time,err_psnr]=test_method1(img,P,run_time,err_psnr,i,j,k,1,flag);
                    % solve problem via TR-ALS
                    if k==1
                        [run_time,err_psnr]=test_method1(img,P,run_time,err_psnr,i,j,k,2,flag);
                    end
                    % solve problem via SiLRTC-TT
                    [run_time,err_psnr]=test_method1(img,P,run_time,err_psnr,i,j,k,3,flag);
                    % solve problem via LRTC-TNN
                    [run_time,err_psnr]=test_method1(img,P,run_time,err_psnr,i,j,k,4,flag);
                    % solve problem via FBCP
                    [run_time,err_psnr]=test_method1(img,P,run_time,err_psnr,i,j,k,5,flag);
                    % solve problem via HaLRTC
                    [run_time,err_psnr]=test_method1(img,P,run_time,err_psnr,i,j,k,6,flag);
                    % solve problem via STTC
                    [run_time,err_psnr]=test_method1(img,P,run_time,err_psnr,i,j,k,7,flag);
                    % solve problem via TRNNM
                    [run_time,err_psnr]=test_method1(img,P,run_time,err_psnr,i,j,k,8,flag);
                end
            end
        case 6
            img=imread('house.bmp');
            img=double(img);
            for j=1:L
                j
                for k=1:num_trial
                    % sampling
                    P=sampling_uniform(img,sr(j));
                    % solve problem via TR-ADMM
                    [run_time,err_psnr]=test_method1(img,P,run_time,err_psnr,i,j,k,1,flag);
                    % solve problem via TR-ALS
                    if k==1
                        [run_time,err_psnr]=test_method1(img,P,run_time,err_psnr,i,j,k,2,flag);
                    end
                    % solve problem via SiLRTC-TT
                    [run_time,err_psnr]=test_method1(img,P,run_time,err_psnr,i,j,k,3,flag);
                    % solve problem via LRTC-TNN
                    [run_time,err_psnr]=test_method1(img,P,run_time,err_psnr,i,j,k,4,flag);
                    % solve problem via FBCP
                    [run_time,err_psnr]=test_method1(img,P,run_time,err_psnr,i,j,k,5,flag);
                    % solve problem via HaLRTC
                    [run_time,err_psnr]=test_method1(img,P,run_time,err_psnr,i,j,k,6,flag);
                    % solve problem via STTC
                    [run_time,err_psnr]=test_method1(img,P,run_time,err_psnr,i,j,k,7,flag);
                    % solve problem via TRNNM
                    [run_time,err_psnr]=test_method1(img,P,run_time,err_psnr,i,j,k,8,flag);
                end
            end
        case 7
            img=imread('airplane.bmp');
            img=double(img);
            for j=1:L
                j
                for k=1:num_trial
                    % sampling
                    P=sampling_uniform(img,sr(j));
                    % solve problem via TR-ADMM
                    [run_time,err_psnr]=test_method1(img,P,run_time,err_psnr,i,j,k,1,flag);
                    % solve problem via TR-ALS
                    if k==1
                        [run_time,err_psnr]=test_method1(img,P,run_time,err_psnr,i,j,k,2,flag);
                    end
                    % solve problem via SiLRTC-TT
                    [run_time,err_psnr]=test_method1(img,P,run_time,err_psnr,i,j,k,3,flag);
                    % solve problem via LRTC-TNN
                    [run_time,err_psnr]=test_method1(img,P,run_time,err_psnr,i,j,k,4,flag);
                    % solve problem via FBCP
                    [run_time,err_psnr]=test_method1(img,P,run_time,err_psnr,i,j,k,5,flag);
                    % solve problem via HaLRTC
                    [run_time,err_psnr]=test_method1(img,P,run_time,err_psnr,i,j,k,6,flag);
                    % solve problem via STTC
                    [run_time,err_psnr]=test_method1(img,P,run_time,err_psnr,i,j,k,7,flag);
                    % solve problem via TRNNM
                    [run_time,err_psnr]=test_method1(img,P,run_time,err_psnr,i,j,k,8,flag);
                end
            end
        case 8
            load Einstein
            img=Data;
            img=double(img);
            for j=1:L
                j
                for k=1:num_trial
                    % sampling
                    P=sampling_uniform(img,sr(j));
                    % solve problem via TR-ADMM
                    [run_time,err_psnr]=test_method1(img,P,run_time,err_psnr,i,j,k,1,flag);
                    % solve problem via TR-ALS
                    if k==1
                        [run_time,err_psnr]=test_method1(img,P,run_time,err_psnr,i,j,k,2,flag);
                    end
                    % solve problem via SiLRTC-TT
                    [run_time,err_psnr]=test_method1(img,P,run_time,err_psnr,i,j,k,3,flag);
                    % solve problem via LRTC-TNN
                    [run_time,err_psnr]=test_method1(img,P,run_time,err_psnr,i,j,k,4,flag);
                    % solve problem via FBCP
                    [run_time,err_psnr]=test_method1(img,P,run_time,err_psnr,i,j,k,5,flag);
                    % solve problem via HaLRTC
                    [run_time,err_psnr]=test_method1(img,P,run_time,err_psnr,i,j,k,6,flag);
                    % solve problem via STTC
                    [run_time,err_psnr]=test_method1(img,P,run_time,err_psnr,i,j,k,7,flag);
                    % solve problem via TRNNM
                    [run_time,err_psnr]=test_method1(img,P,run_time,err_psnr,i,j,k,8,flag);
                end
            end
    end
end
warning on
%%
if flag
    save run_time1_VDT run_time
    save err_psnr1_VDT err_psnr
else
    save run_time1 run_time
    save err_psnr1 err_psnr
end
if num_trial==1
    run_time_s=cellfun(@squeeze,run_time,'UniformOutput',false);
    err_psnr_s=cellfun(@squeeze,err_psnr,'UniformOutput',false);
else
    run_time_s=cellfun(@mean,run_time,'UniformOutput',false);
    run_time_s=cellfun(@squeeze,run_time_s,'UniformOutput',false);
    err_psnr_s=cellfun(@mean,err_psnr,'UniformOutput',false);
    err_psnr_s=cellfun(@squeeze,err_psnr_s,'UniformOutput',false);
end
if flag
    save run_time_s1_VDT run_time_s
    save err_psnr_s1_VDT err_psnr_s
else
    save run_time_s1 run_time_s
    save err_psnr_s1 err_psnr_s
end
%%
run_time_c=cat(3,run_time_s{1},run_time_s{2},run_time_s{3},run_time_s{4},...
    run_time_s{5},run_time_s{6},run_time_s{7},run_time_s{8});
err_psnr_c=cat(3,err_psnr_s{1},err_psnr_s{2},err_psnr_s{3},err_psnr_s{4},...
    err_psnr_s{5},err_psnr_s{6},err_psnr_s{7},err_psnr_s{8});
%% visualizaion
str=["kodim04","peppers","sailboat","lena","barbara","house","airplane","Einstein"];
figure;
for i=1:2
    for j=1:4
        subplot(2,4,4*(i-1)+j);
        p=plot(sr,squeeze(err_psnr_c(4*(i-1)+j,:,:)));
        p(1).Marker='p';
        p(2).Marker='o';
        p(3).Marker='*';
        p(4).Marker='x';
        p(5).Marker='s';
        p(6).Marker='d';
        p(7).Marker='^';
        p(8).Marker='+';
        xlabel('SR');
        ylabel('PSNR (dB)');
        title(str(4*(i-1)+j));
    end
end
hl=legend({'TRBU','TR-ALS','SiLRTC-TT','LRTC-TNN',...
    'FBCP','HaLRTC','STTC','TRNNM'});
set(hl,'box','off');
figure;
for i=1:2
    for j=1:4
        subplot(2,4,4*(i-1)+j);
        p=semilogy(sr,squeeze(run_time_c(4*(i-1)+j,:,:)));
        p(1).Marker='p';
        p(2).Marker='o';
        p(3).Marker='*';
        p(4).Marker='x';
        p(5).Marker='s';
        p(6).Marker='d';
        p(7).Marker='^';
        p(8).Marker='+';
        xlabel('SR');
        ylabel('CPUtime (s)');
        title(str(4*(i-1)+j));
    end
end