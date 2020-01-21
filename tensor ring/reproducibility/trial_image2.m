clc
clear
close all
%% input arguments
num_trial=1;
flag=true;
%% image tests
X=cell(2,8);
err_psnr=cell(1,8);
run_time=cell(1,8);
% load X2_VDT
% load err_psnr2_VDT
% load run_time2_VDT
warning off
for i=1:2
    i
    switch i
        case 1
            img=imread('house.bmp');
            img=double(img);
            I=imread('printedtext.png');
            BW=imbinarize(I,'adaptive','ForegroundPolarity','dark','Sensitivity',0.4);
            P=imresize(BW,[256,256]);
            P=cat(3,P,P,P);
            P=double(P);
            for k=1:num_trial
                %solve problem via TR-ADMM
                [X,run_time,err_psnr]=test_method2(img,P,X,run_time,err_psnr,i,k,1,flag);
                % solve problem via TR-ALS
                [X,run_time,err_psnr]=test_method2(img,P,X,run_time,err_psnr,i,k,2,flag);
                % solve problem via SiLRTC-TT
                [X,run_time,err_psnr]=test_method2(img,P,X,run_time,err_psnr,i,k,3,flag);
                % solve problem via FP-LRTC
                [X,run_time,err_psnr]=test_method2(img,P,X,run_time,err_psnr,i,k,4,flag);
                % solve problem via FBCP
                [X,run_time,err_psnr]=test_method2(img,P,X,run_time,err_psnr,i,k,5,flag);
                % solve problem via HaLRTC
                [X,run_time,err_psnr]=test_method2(img,P,X,run_time,err_psnr,i,k,6,flag);
                % solve problem via STTC
                [X,run_time,err_psnr]=test_method2(img,P,X,run_time,err_psnr,i,k,7,flag);
                % solve problem via TRNNM
                [X,run_time,err_psnr]=test_method2(img,P,X,run_time,err_psnr,i,k,8,flag);
            end
        case 2
            img1=imread('llama.jpg');
            img=zeros([256,384,3]);
            for n=1:3
                img(:,:,n)=imresize(img1(:,:,n),[256,384]);
            end
            img=double(img);
            load P
            for k=1:num_trial
                %solve problem via TR-ADMM
                [X,run_time,err_psnr]=test_method2(img,P,X,run_time,err_psnr,i,k,1,flag);
                % solve problem via TR-ALS
                [X,run_time,err_psnr]=test_method2(img,P,X,run_time,err_psnr,i,k,2,flag);
                % solve problem via SiLRTC-TT
                [X,run_time,err_psnr]=test_method2(img,P,X,run_time,err_psnr,i,k,3,flag);
                % solve problem via FP-LRTC
                [X,run_time,err_psnr]=test_method2(img,P,X,run_time,err_psnr,i,k,4,flag);
                % solve problem via FBCP
                [X,run_time,err_psnr]=test_method2(img,P,X,run_time,err_psnr,i,k,5,flag);
                % solve problem via HaLRTC
                [X,run_time,err_psnr]=test_method2(img,P,X,run_time,err_psnr,i,k,6,flag);
                % solve problem via STTC
                [X,run_time,err_psnr]=test_method2(img,P,X,run_time,err_psnr,i,k,7,flag);
                % solve problem via TRNNM
                [X,run_time,err_psnr]=test_method2(img,P,X,run_time,err_psnr,i,k,8,flag);
            end
    end
end
warning on
%%
if flag
    save run_time2_VDT run_time
    save err_psnr2_VDT err_psnr
    save X2_VDT X
else
    save run_time2 run_time
    save err_psnr2 err_psnr
    save X2 X
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
img=imread('house.bmp');
img=double(img);
I=imread('printedtext.png');
BW=imbinarize(I,'adaptive','ForegroundPolarity','dark','Sensitivity',0.4);
P=imresize(BW,[256,256]);
P=cat(3,P,P,P);
P=double(P);
figure;
imshow(uint8(img),'Border','tight');
saveas(gcf,'original_house.png');
figure;
imshow(uint8(P.*img),'Border','tight');
saveas(gcf,'observed_house.png');
for i=1:8
    figure;
    imshow(uint8(X{1,i}),'Border','tight');
    if flag
        saveas(gcf,['house_' str{i} '_VDT.png']);
    else
        saveas(gcf,['house_' str{i} '.png']);
    end
end
img1=imread('llama.jpg');
img=zeros([256,384,3]);
for l=1:3
    img(:,:,l)=imresize(img1(:,:,l),[256,384]);
end
siz=size(img);
img=double(img);
load P
figure;
imshow(uint8(img),'Border','tight');
saveas(gcf,'original_llama.png');
figure;
imshow(uint8(P.*img),'Border','tight');
saveas(gcf,'observed_llama.png');
for i=1:8
    figure;
    imshow(uint8(X{2,i}),'Border','tight');
    if flag
        saveas(gcf,['llama_' str{i} '_VDT.png']);
    else
        saveas(gcf,['llama_' str{i} '.png']);
    end
end