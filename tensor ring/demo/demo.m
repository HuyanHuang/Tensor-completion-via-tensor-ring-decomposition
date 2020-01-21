clc
clear
close all
%% input arguments
TRr=7;
sr=0.1;
img=imread('lena.bmp');
%% sampling
siz=size(img);
img=double(img);
P=sampling_uniform(img,sr);
%% ket augmentation--covert original(low-order) tensor to high-order tensor
I=[2*ones(1,16) 3];
order=[1 9 2 10 3 11 4 12 5 13 6 14 7 15 8 16 17];
J=[4*ones(1,8) 3];
tnsr=l2h(img,I,order,J);
P=l2h(P,I,order,J);
%% solve problem via ADMM/ALS
[x,RC,RE,run_time]=TR_ADMM(tnsr,P,10^-3.7,true);
% [x,RC,RE,run_time]=TR_ALS(tnsr,P,TRr,true);
%% evaluation
x=h2l(x,I,order,siz);
P=h2l(P,I,order,siz);
err_re=RE(end)
err_psnr=psnr(uint8(x),uint8(img))
% err_ssim=ssim(uint8(reshape(x,siz)),uint8(img))
%% deaugmentation--covert high-order tensor to original(low-order) tensor
img_recovered=uint8(x);
img_observed=uint8(P.*img);
%% visualize the results
figure;
imshow(img_observed,'border','tight');
figure;
imshow(img_recovered,'border','tight');