function [run_time,err_psnr]=test_method1(img,P,run_time,err_psnr,i,j,k,method,flag)
siz=size(img);
switch method
    case 1
        %% solve problem via TR-ADMM (high)
        if flag
            if i==1
                I=[2*ones(1,8),3,2*ones(1,9),3];
                order=[1:(length(I)-1)/2;(length(I)+1)/2:length(I)-1];
                order=[order(:);length(I)]';
                J=[4*ones(1,8),6,3];
                tnsr=l2h(img,I,order,J);
                P=l2h(P,I,order,J);
            elseif i==8
                I=[2,2,2,3,5,5,5,5,3,2,2,2,3];
                order=[1:(length(I)-1)/2;(length(I)+1)/2:length(I)-1];
                order=[order(:);length(I)]';
                J=[10,10,6,6,10,10,3];
                tnsr=l2h(img,I,order,J);
                P=l2h(P,I,order,J);
            else
                I=[2*ones(1,16),3];
                order=[1:(length(I)-1)/2;(length(I)+1)/2:length(I)-1];
                order=[order(:);length(I)]';
                J=[4*ones(1,8),3];
                tnsr=l2h(img,I,order,J);
                P=l2h(P,I,order,J);
            end
            [x,~,~,run_time{method}(k,i,j)]=TR_ADMM(tnsr,P,10^-3.7,false);
            x=h2l(x,I,order,siz);
        else
            [x,~,~,run_time{method}(k,i,j)]=TR_ADMM(img,P,10^-3.7,false);
        end
        err_psnr{method}(k,i,j)=psnr(uint8(x),uint8(img));
    case 2
        %% solve problem via TR-ALS
        if flag
            if i==1
                I=[2*ones(1,8),3,2*ones(1,9),3];
                order=[1:(length(I)-1)/2;(length(I)+1)/2:length(I)-1];
                order=[order(:);length(I)]';
                J=[4*ones(1,8),6,3];
                tnsr=l2h(img,I,order,J);
                P=l2h(P,I,order,J);
                TRr=12;
            elseif i==8
                I=[2,2,2,3,5,5,5,5,3,2,2,2,3];
                order=[1:(length(I)-1)/2;(length(I)+1)/2:length(I)-1];
                order=[order(:);length(I)]';
                J=[10,10,6,6,10,10,3];
                tnsr=l2h(img,I,order,J);
                P=l2h(P,I,order,J);
                TRr=14;
            else
                I=[2*ones(1,16),3];
                order=[1:(length(I)-1)/2;(length(I)+1)/2:length(I)-1];
                order=[order(:);length(I)]';
                J=[4*ones(1,8),3];
                tnsr=l2h(img,I,order,J);
                P=l2h(P,I,order,J);
                TRr=14;
            end
            [x,~,~,run_time{method}(k,i,j)]=TR_ALS(tnsr,P,TRr,false);
            x=h2l(x,I,order,siz);
        else
            if i==1
                TRr=12;
            else
                TRr=14;
            end
            [x,~,~,run_time{method}(k,i,j)]=TR_ALS(img,P,TRr,false);
        end
        err_psnr{method}(k,i,j)=psnr(uint8(x),uint8(img));
    case 3
        %% solve problem via SiLRTC-TT
        if flag
            if i==1
                I=[2*ones(1,8),3,2*ones(1,9),3];
                order=[1:(length(I)-1)/2;(length(I)+1)/2:length(I)-1];
                order=[order(:);length(I)]';
                J=[4*ones(1,8),6,3];
                tnsr=l2h(img,I,order,J);
                P=l2h(P,I,order,J);
            elseif i==8
                I=[2,2,2,3,5,5,5,5,3,2,2,2,3];
                order=[1:(length(I)-1)/2;(length(I)+1)/2:length(I)-1];
                order=[order(:);length(I)]';
                J=[10,10,6,6,10,10,3];
                tnsr=l2h(img,I,order,J);
                P=l2h(P,I,order,J);
            else
                tnsr=CastImageAsKet(img,[4*ones(1,8),3],2,2);
                P=CastImageAsKet(P,[4*ones(1,8),3],2,2);
            end
            t=cputime;
            [~,ranktube]=SVD_MPS_Rank_Estimation(P.*tnsr,1e-2);
            x=SiLRTC_TT(tnsr,find(P),ranktube,0.01*ranktube,500,1e-5);
            run_time{method}(k,i,j)=cputime-t;
            if i~=1 && i~=8
                x=CastKet2Image(x,256,256,2,2);
            else
                x=h2l(x,I,order,siz);
            end
        else
            t=cputime;
            [~,ranktube]=SVD_MPS_Rank_Estimation(P.*img,1e-2);
            x=SiLRTC_TT(img,find(P),ranktube,0.05*ranktube,500,1e-5);
            run_time{method}(k,i,j)=cputime-t;
        end
        err_psnr{method}(k,i,j)=psnr(uint8(x),uint8(img));
    case 4
        %% solve problem via LRTC-TNN
        if flag
            if i==1
                I=[2*ones(1,8),3,2*ones(1,9),3];
                order=[1:(length(I)-1)/2;(length(I)+1)/2:length(I)-1];
                order=[order(:);length(I)]';
                J=[4*ones(1,8),6,3];
                tnsr=l2h(img,I,order,J);
                P=l2h(P,I,order,J);
            elseif i==8
                I=[2,2,2,3,5,5,5,5,3,2,2,2,3];
                order=[1:(length(I)-1)/2;(length(I)+1)/2:length(I)-1];
                order=[order(:);length(I)]';
                J=[10,10,6,6,10,10,3];
                tnsr=l2h(img,I,order,J);
                P=l2h(P,I,order,J);
            else
                I=[2*ones(1,16),3];
                order=[1:(length(I)-1)/2;(length(I)+1)/2:length(I)-1];
                order=[order(:);length(I)]';
                J=[4*ones(1,8),3];
                tnsr=l2h(img,I,order,J);
                P=l2h(P,I,order,J);
            end
            tnsr=reshape(tnsr,siz);
            P=reshape(P,siz);
            t=cputime;
            x=lrtc_tnn(P.*tnsr,find(P),[]);
            run_time{method}(k,i,j)=cputime-t;
        else
            t=cputime;
            x=lrtc_tnn(P.*img,find(P),[]);
            run_time{method}(k,i,j)=cputime-t;
        end
        err_psnr{method}(k,i,j)=psnr(uint8(x),uint8(img));
    case 5
        %% solve problem via FBCP
        if flag
            if i==1
                I=[2*ones(1,8),3,2*ones(1,9),3];
                order=[1:(length(I)-1)/2;(length(I)+1)/2:length(I)-1];
                order=[order(:);length(I)]';
                J=[4*ones(1,8),6,3];
                tnsr=l2h(img,I,order,J);
                P=l2h(P,I,order,J);
                CPr=50;
            elseif i==8
                I=[2,2,2,3,5,5,5,5,3,2,2,2,3];
                order=[1:(length(I)-1)/2;(length(I)+1)/2:length(I)-1];
                order=[order(:);length(I)]';
                J=[10,10,6,6,10,10,3];
                tnsr=l2h(img,I,order,J);
                P=l2h(P,I,order,J);
                CPr=60;
            else
                I=[2*ones(1,16),3];
                order=[1:(length(I)-1)/2;(length(I)+1)/2:length(I)-1];
                order=[order(:);length(I)]';
                J=[4*ones(1,8),3];
                tnsr=l2h(img,I,order,J);
                P=l2h(P,I,order,J);
                CPr=100;
            end
            t=cputime;
            model=BCPF_IC(P.*tnsr,'obs',P,'init','rand','maxRank',CPr,'maxiters',100,...
                'tol',1e-4,'dimRed',1,'verbose',0);
            run_time{method}(k,i,j)=cputime-t;
            x=double(model.X);
        else
            if i==1
                CPr=50;
            elseif i==8
                CPr=60;
            else
                CPr=100;
            end
            t=cputime;
            model=BCPF_IC(P.*img,'obs',P,'init','rand','maxRank',CPr,'maxiters',100,...
                'tol',1e-4,'dimRed',1,'verbose',0);
            run_time{method}(k,i,j)=cputime-t;
            x=double(model.X);
        end
        err_psnr{method}(k,i,j)=psnr(uint8(x),uint8(img));
    case 6
        %% solve problem via HaLRTC
        if flag
            if i==1
                I=[2*ones(1,8),3,2*ones(1,9),3];
                order=[1:(length(I)-1)/2;(length(I)+1)/2:length(I)-1];
                order=[order(:);length(I)]';
                J=[4*ones(1,8),6,3];
                tnsr=l2h(img,I,order,J);
                P=l2h(P,I,order,J);
            elseif i==8
                I=[2,2,2,3,5,5,5,5,3,2,2,2,3];
                order=[1:(length(I)-1)/2;(length(I)+1)/2:length(I)-1];
                order=[order(:);length(I)]';
                J=[10,10,6,6,10,10,3];
                tnsr=l2h(img,I,order,J);
                P=l2h(P,I,order,J);
            else
                I=[2*ones(1,16),3];
                order=[1:(length(I)-1)/2;(length(I)+1)/2:length(I)-1];
                order=[order(:);length(I)]';
                J=[4*ones(1,8),3];
                tnsr=l2h(img,I,order,J);
                P=l2h(P,I,order,J);
            end
            t=cputime;
            x=HaLRTC(tnsr,logical(P),ones(1,ndims(tnsr))/ndims(tnsr),1e-6,500,1e-5,P.*img);
            run_time{method}(k,i,j)=cputime-t;
        else
            t=cputime;
            x=HaLRTC(img,logical(P),ones(1,length(siz))/length(siz),1e-6,500,1e-5,P.*img);
            run_time{method}(k,i,j)=cputime-t;
        end
        err_psnr{method}(k,i,j)=psnr(uint8(x),uint8(img));
    case 7
        %% solve problem via STTC
        if flag
            if i==1
                I=[2*ones(1,8),3,2*ones(1,9),3];
                order=[1:(length(I)-1)/2;(length(I)+1)/2:length(I)-1];
                order=[order(:);length(I)]';
                J=[4*ones(1,8),6,3];
                tnsr=l2h(img,I,order,J);
                P=l2h(P,I,order,J);
            elseif i==8
                I=[2,2,2,3,5,5,5,5,3,2,2,2,3];
                order=[1:(length(I)-1)/2;(length(I)+1)/2:length(I)-1];
                order=[order(:);length(I)]';
                J=[10,10,6,6,10,10,3];
                tnsr=l2h(img,I,order,J);
                P=l2h(P,I,order,J);
            else
                I=[2*ones(1,16),3];
                order=[1:(length(I)-1)/2;(length(I)+1)/2:length(I)-1];
                order=[order(:);length(I)]';
                J=[4*ones(1,8),3];
                tnsr=l2h(img,I,order,J);
                P=l2h(P,I,order,J);
            end
            tnsr=reshape(tnsr,siz);
            P=reshape(P,siz);
            t=cputime;
            x=Smoothlowrank_TV12(P.*tnsr,find(P),0,P.*tnsr,[10,10,0]);
            run_time{method}(k,i,j)=cputime-t;
        else
            t=cputime;
            x=Smoothlowrank_TV12(P.*img,find(P),0,P.*img,[10,10,0]);
            run_time{method}(k,i,j)=cputime-t;
        end
        err_psnr{method}(k,i,j)=psnr(uint8(x),uint8(img));
    case 8
        %% solve problem via TRNNM
        if flag
            if i==1
                I=[2*ones(1,8),3,2*ones(1,9),3];
                order=[1:(length(I)-1)/2;(length(I)+1)/2:length(I)-1];
                order=[order(:);length(I)]';
                J=[4*ones(1,8),6,3];
                tnsr=l2h(img,I,order,J);
                P=l2h(P,I,order,J);
            elseif i==8
                I=[2,2,2,3,5,5,5,5,3,2,2,2,3];
                order=[1:(length(I)-1)/2;(length(I)+1)/2:length(I)-1];
                order=[order(:);length(I)]';
                J=[10,10,6,6,10,10,3];
                tnsr=l2h(img,I,order,J);
                P=l2h(P,I,order,J);
            else
                I=[2*ones(1,16),3];
                order=[1:(length(I)-1)/2;(length(I)+1)/2:length(I)-1];
                order=[order(:);length(I)]';
                J=[4*ones(1,8),3];
                tnsr=l2h(img,I,order,J);
                P=l2h(P,I,order,J);
            end
            [x,~,~,run_time{method}(k,i,j)]=TRNNM(tnsr,P,10^-5,false);
            x=h2l(x,I,order,siz);
        else
            [x,~,~,run_time{method}(k,i,j)]=TRNNM(img,P,10^-5,false);
        end
        err_psnr{method}(k,i,j)=psnr(uint8(x),uint8(img));
end
end