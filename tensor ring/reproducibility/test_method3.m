function [X,run_time,err_psnr]=test_method3(vid,P,X,run_time,err_psnr,i,k,method,flag)
siz=size(vid);
switch method
    case 1
        %% solve problem via TR-ADMM
        if flag
            if i==1
                I=[2,2,2,2,5,4,4,2,2,2,3,4,5,5];
                order=[1:5;6:10];
                order=[order(:);(11:length(I))']';
                J=[8,8,4,4,10,3,4,5,5];
                tnsr=l2h(vid,I,order,J);
                P=l2h(P,I,order,J);
            else
                I=[2,2,2,3,3,4,4,2,2,2,3,4,5,5];
                order=[1:5;6:10];
                order=[order(:);(11:length(I))']';
                J=[8,8,4,6,6,3,4,5,5];
                tnsr=l2h(vid,I,order,J);
                P=l2h(P,I,order,J);
            end
            [x,~,~,run_time{method}(k,i)]=TR_ADMM(tnsr,P,1e-3,false);
            x=h2l(x,I,order,siz);
        else
            [x,~,~,run_time{method}(k,i)]=TR_ADMM(vid,P,1e-3,false);
        end
        err_psnr{method}(k,i)=psnr(uint8(x),uint8(vid));
        X{i,method}=x;
    case 2
        %% solve problem via TR-ALS
        if flag
            if i==1
                I=[2,2,2,2,5,4,4,2,2,2,3,4,5,5];
                order=[1:5;6:10];
                order=[order(:);(11:length(I))']';
                J=[8,8,4,4,10,3,4,5,5];
                tnsr=l2h(vid,I,order,J);
                P=l2h(P,I,order,J);
            else
                I=[2,2,2,3,3,4,4,2,2,2,3,4,5,5];
                order=[1:5;6:10];
                order=[order(:);(11:length(I))']';
                J=[8,8,4,6,6,3,4,5,5];
                tnsr=l2h(vid,I,order,J);
                P=l2h(P,I,order,J);
            end
            TRr=12;
            [x,~,~,run_time{method}(k,i)]=TR_ALS(tnsr,P,TRr,false);
            x=h2l(x,I,order,siz);
        else
            TRr=18;
            [x,~,~,run_time{method}(k,i)]=TR_ALS(vid,P,TRr,false);
        end
        err_psnr{method}(k,i)=psnr(uint8(x),uint8(vid));
        X{i,method}=x;
    case 3
        %% solve problem via SiLRTC-TT
        if flag
            if i==1
                I=[2,2,2,2,5,4,4,2,2,2,3,4,5,5];
                order=[1:5;6:10];
                order=[order(:);(11:length(I))']';
                J=[8,8,4,4,10,3,4,5,5];
                tnsr=l2h(vid,I,order,J);
                P=l2h(P,I,order,J);
            else
                I=[2,2,2,3,3,4,4,2,2,2,3,4,5,5];
                order=[1:5;6:10];
                order=[order(:);(11:length(I))']';
                J=[8,8,4,6,6,3,4,5,5];
                tnsr=l2h(vid,I,order,J);
                P=l2h(P,I,order,J);
            end
            t=cputime;
            [~,ranktube]=SVD_MPS_Rank_Estimation(P.*tnsr,1e-2);
            x=SiLRTC_TT(tnsr,find(P),ranktube,0.01*ranktube,500,1e-5);
            run_time{method}(k,i)=cputime-t;
            x=h2l(x,I,order,siz);
        else
            t=cputime;
            [~,ranktube]=SVD_MPS_Rank_Estimation(P.*vid,1e-2);
            x=SiLRTC_TT(vid,find(P),ranktube,0.05*ranktube,500,1e-5);
            run_time{method}(k,i)=cputime-t;
        end
        err_psnr{method}(k,i)=psnr(uint8(x),uint8(vid));
        X{i,method}=x;
    case 4
        %% solve problem via LRTC-TNN
        if flag
            if i==1
                I=[2,2,2,2,5,4,4,2,2,2,3,4,5,5];
                order=[1:5;6:10];
                order=[order(:);(11:length(I))']';
                J=[8,8,4,4,10,3,4,5,5];
                tnsr=l2h(vid,I,order,J);
                P=l2h(P,I,order,J);
            else
                I=[2,2,2,3,3,4,4,2,2,2,3,4,5,5];
                order=[1:5;6:10];
                order=[order(:);(11:length(I))']';
                J=[8,8,4,6,6,3,4,5,5];
                tnsr=l2h(vid,I,order,J);
                P=l2h(P,I,order,J);
            end
            tnsr=reshape(tnsr,[siz(1),siz(2),prod(siz(3:end))]);
            P=reshape(P,[siz(1),siz(2),prod(siz(3:end))]);
            t=cputime;
            x=lrtc_tnn(P.*tnsr,find(P),[]);
            run_time{method}(k,i)=cputime-t;
        else
            vid=reshape(vid,[siz(1),siz(2),prod(siz(3:end))]);
            P=reshape(P,[siz(1),siz(2),prod(siz(3:end))]);
            t=cputime;
            x=lrtc_tnn(P.*vid,find(P),[]);
            run_time{method}(k,i)=cputime-t;
        end
        x=reshape(x,siz);
        vid=reshape(vid,siz);
        err_psnr{method}(k,i)=psnr(uint8(x),uint8(vid));
        X{i,method}=x;
    case 5
        %% solve problem via FBCP
        if flag
            if i==1
                I=[2,2,2,2,5,4,4,2,2,2,3,4,5,5];
                order=[1:5;6:10];
                order=[order(:);(11:length(I))']';
                J=[8,8,4,4,10,3,4,5,5];
                tnsr=l2h(vid,I,order,J);
                P=l2h(P,I,order,J);
            else
                I=[2,2,2,3,3,4,4,2,2,2,3,4,5,5];
                order=[1:5;6:10];
                order=[order(:);(11:length(I))']';
                J=[8,8,4,6,6,3,4,5,5];
                tnsr=l2h(vid,I,order,J);
                P=l2h(P,I,order,J);
            end
            CPr=40;
            t=cputime;
            model=BCPF_IC(P.*tnsr,'obs',P,'init','rand','maxRank',CPr,'maxiters',100,...
                'tol',1e-4,'dimRed',1,'verbose',0);
            run_time{method}(k,i)=cputime-t;
            x=double(model.X);
        else
            CPr=40;
            t=cputime;
            model=BCPF_IC(P.*vid,'obs',P,'init','rand','maxRank',CPr,'maxiters',100,...
                'tol',1e-4,'dimRed',1,'verbose',0);
            run_time{method}(k,i)=cputime-t;
            x=double(model.X);
        end
        err_psnr{method}(k,i)=psnr(uint8(x),uint8(vid));
        X{i,method}=x;
    case 6
        %% solve problem via HaLRTC
        t=cputime;
        x=HaLRTC(vid,logical(P),ones(1,length(siz))/length(siz),1e-6,500,1e-5,P.*vid);
        run_time{method}(k,i)=cputime-t;
        err_psnr{method}(k,i)=psnr(uint8(x),uint8(vid));
        X{i,method}=x;
    case 7
        %% solve problem via STTC
        if flag
            if i==1
                I=[2,2,2,2,5,4,4,2,2,2,3,4,5,5];
                order=[1:5;6:10];
                order=[order(:);(11:length(I))']';
                J=[8,8,4,4,10,3,4,5,5];
                tnsr=l2h(vid,I,order,J);
                P=l2h(P,I,order,J);
            else
                I=[2,2,2,3,3,4,4,2,2,2,3,4,5,5];
                order=[1:5;6:10];
                order=[order(:);(11:length(I))']';
                J=[8,8,4,6,6,3,4,5,5];
                tnsr=l2h(vid,I,order,J);
                P=l2h(P,I,order,J);
            end
            tnsr=reshape(tnsr,siz);
            P=reshape(P,siz);
            t=cputime;
            x=Smoothlowrank_TV_video21(P.*tnsr,find(P),0,P.*tnsr,[20,20,0,20]);
            run_time{method}(k,i)=cputime-t;
        else
            t=cputime;
            x=Smoothlowrank_TV_video21(P.*vid,find(P),0,P.*vid,[20,20,0,20]);
            run_time{method}(k,i)=cputime-t;
        end
        err_psnr{method}(k,i)=psnr(uint8(x),uint8(vid));
        X{i,method}=x;
    case 8
        %% solve problem via TRNNM
        if flag
            if i==1
                I=[2,2,2,2,5,4,4,2,2,2,3,4,5,5];
                order=[1:5;6:10];
                order=[order(:);(11:length(I))']';
                J=[8,8,4,4,10,3,4,5,5];
                tnsr=l2h(vid,I,order,J);
                P=l2h(P,I,order,J);
            else
                I=[2,2,2,3,3,4,4,2,2,2,3,4,5,5];
                order=[1:5;6:10];
                order=[order(:);(11:length(I))']';
                J=[8,8,4,6,6,3,4,5,5];
                tnsr=l2h(vid,I,order,J);
                P=l2h(P,I,order,J);
            end
            [x,~,~,run_time{method}(k,i)]=TRNNM(tnsr,P,1e-5,false);
            x=h2l(x,I,order,siz);
        else
            [x,~,~,run_time{method}(k,i)]=TRNNM(vid,P,1e-5,false);
        end
        err_psnr{method}(k,i)=psnr(uint8(x),uint8(vid));
        X{i,method}=x;
end
end