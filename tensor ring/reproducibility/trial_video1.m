clc
clear
close all
%% input arguments
sr=0.1;
num_trial=1;
flag=true;
%% image tests
X=cell(2,8);
err_psnr=cell(1,8);
run_time=cell(1,8);
% load X3_VDT
% load err_psnr3_VDT
% load run_time3_VDT
warning off
for i=1:2
    i
    switch i
        case 1
            load vid_explosion
            vid1=double(vid);
            P1=sampling_uniform(vid1,sr);
            for k=1:num_trial
                % solve problem via TR-ADMM
                 [X,run_time,err_psnr]=test_method3(vid1,P1,X,run_time,err_psnr,i,k,1,flag);
                % solve problem via TR-ALS
                [X,run_time,err_psnr]=test_method3(vid1,P1,X,run_time,err_psnr,i,k,2,flag);
                % solve problem via SiLRTC-TT
                [X,run_time,err_psnr]=test_method3(vid1,P1,X,run_time,err_psnr,i,k,3,flag);
                % solve problem via LRTC-TNN
                [X,run_time,err_psnr]=test_method3(vid1,P1,X,run_time,err_psnr,i,k,4,flag);
                % solve problem via FBCP
                [X,run_time,err_psnr]=test_method3(vid1,P1,X,run_time,err_psnr,i,k,5,flag);
                % solve problem via HaLRTC
                [X,run_time,err_psnr]=test_method3(vid1,P1,X,run_time,err_psnr,i,k,6,flag);
                % solve problem via STTC
                [X,run_time,err_psnr]=test_method3(vid1,P1,X,run_time,err_psnr,i,k,7,flag);
                % solve problem via TRNNM
                [X,run_time,err_psnr]=test_method3(vid1,P1,X,run_time,err_psnr,i,k,8,flag);
            end
        case 2
            load vid_cock
            vid2=double(vid);
            P2=sampling_uniform(vid2,sr);
            for k=1:num_trial
                % solve problem via TR-ADMM
                 [X,run_time,err_psnr]=test_method3(vid2,P2,X,run_time,err_psnr,i,k,1,flag);
                % solve problem via TR-ALS
                [X,run_time,err_psnr]=test_method3(vid2,P2,X,run_time,err_psnr,i,k,2,flag);
                % solve problem via SiLRTC-TT
                [X,run_time,err_psnr]=test_method3(vid2,P2,X,run_time,err_psnr,i,k,3,flag);
                % solve problem via LRTC-TNN
                [X,run_time,err_psnr]=test_method3(vid2,P2,X,run_time,err_psnr,i,k,4,flag);
                % solve problem via FBCP
                [X,run_time,err_psnr]=test_method3(vid2,P2,X,run_time,err_psnr,i,k,5,flag);
                % solve problem via HaLRTC
                [X,run_time,err_psnr]=test_method3(vid2,P2,X,run_time,err_psnr,i,k,6,flag);
                % solve problem via STTC
                [X,run_time,err_psnr]=test_method3(vid2,P2,X,run_time,err_psnr,i,k,7,flag);
                % solve problem via TRNNM
                [X,run_time,err_psnr]=test_method3(vid2,P2,X,run_time,err_psnr,i,k,8,flag);
            end
    end
end
warning on
%%
if flag
save run_time3_VDT run_time
save err_psnr3_VDT err_psnr
save X3_VDT X
else
    save run_time3 run_time
save err_psnr3 err_psnr
save X3 X
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
run_time_c=cat(2,run_time_s{1},run_time_s{2},run_time_s{3},run_time_s{4},...
    run_time_s{5},run_time_s{6},run_time_s{7},run_time_s{8});
err_psnr_c=cat(2,err_psnr_s{1},err_psnr_s{2},err_psnr_s{3},err_psnr_s{4},...
    err_psnr_s{5},err_psnr_s{6},err_psnr_s{7},err_psnr_s{8});
%% visualizaion
str={'TRBU','TR-ALS','SiLRTC-TT','LRTC-TNN','FBCP','HaLRTC','STTC','TRNNM'};
% figure;
% imshow(uint8(vid1(:,:,:,end)),'Border','tight');
% saveas(gcf,'original_explosion.png');
% figure;
% imshow(uint8(P1(:,:,:,end).*vid1(:,:,:,end)),'Border','tight');
% saveas(gcf,'observed_explosion.png');
for i=1:8
    figure;
    imshow(uint8(X{1,i}(:,:,:,end)),'Border','tight');
    if flag
        saveas(gcf,['explosion_' str{i} '_VDT.png']);
    else
        saveas(gcf,['explosion_' str{i} '.png']);
    end
end
% figure;
% imshow(uint8(vid2(:,:,:,1)),'Border','tight');
% saveas(gcf,'original_cock.png');
% figure;
% imshow(uint8(P2(:,:,:,1).*vid2(:,:,:,1)),'Border','tight');
% saveas(gcf,'observed_cock.png');
for i=1:8
    figure;
    imshow(uint8(X{2,i}(:,:,:,1)),'Border','tight');
    if flag
        saveas(gcf,['cock_' str{i} '_VDT.png']);
    else
        saveas(gcf,['cock_' str{i} '.png']);
    end
end