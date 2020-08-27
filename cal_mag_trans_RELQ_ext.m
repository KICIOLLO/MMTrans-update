function cal_mag_trans_RELQ_ext(str_open,I0_ang,D0_ang,DirRes,SavFor,ext_para)

    % 函数说明： 总强度磁异常 转换 多个参量 (R E L Q)
    
    % Updated by KICIOLLO (sdsun@hust.edu.cn; shidasun.cug@gmail.com; http://kiciollo.github.io/)
    
    % Based on the Matlab scripts (MMTrans) by Gerovska and Araúzo-Bravo (2006)
    % Gerovska, D., and M. J. Araúzo-Bravo, 2006, Calculation of magnitude magnetic transforms 
    %      with high centricity and low dependence on the magnetization vector direction: Geophysics, 
    %      71, no. 5, 121-130, doi: 10.1190/1.2335516.
    
    % 输入参数： str_open       网格数据filename或者是pathname都可以
    %           I0_ang D0_ang  地磁场（正常场）倾角及偏角
    
    %           DirRes         计算结果存放路径
    %           SavFor         计算结果保存格式 最好是 'grd' 这种形式
    
    %           ext_para
    %               = 1   余弦函数扩边
    %               = 2   多项式函数扩边
    
    % 改进历史：
    %         20140317 未保存记录
    %         20140317 22:17 拟加入张壹所编扩边函数，因自带扩边函数ExtCos只扩左边和下边
    %         20140317 22:30 加入扩边函数完成！！
    %         20140319 22:56 在函数Delta2HaxHayZa中不再进行Ta的计算部分，改为后面由xyz求取
    %                        并且求取的X Y Z T都不再主动删除对应的mat文件，所以每次函数调用
    %                        后需要配合使用  delete('X.mat'); 语句
    %         20140320 21:34 matlab2012a 版本，fileparts函数跟2009版本不同，需要注意
    %         20140321 14:51 在实际使用过程中，发现在含噪声时，多项式扩边出现了特别大的严重错误
    %                        因此将在程序中仍加入 余弦函数扩边
    %         20140321 17:06 余弦扩边完成（并将PER改成了20）
    
    %         20170316 16:42 拟将程序更改为 化极 运算程序；当前先不考虑低纬度化极的情况
    %                  20:37 化极程序已经搞定，这里是打算将根据 deltaT 及正常场方向来计算 NSS

    %         20170827 18:00 在NSS计算的基础上，进行修改，完成计算参量 R E L Q 的工作
    
    
    
    FilSep=filesep;
    Per = 20;
    [Dir,Nam,Ext]=fileparts(str_open);
    
    % Test if the directory of the results exists
    if ~isdir(DirRes)
       disp(['The output directory ' DirRes ' does not exist']);
       disp(['It will be created']);
       mkdir(DirRes);
    end % if ~isdir(DirRes)
    
    % Extension results
    ExtRes=['.' SavFor];
    
    % Read data
    ReaFor=Ext;
    switch ReaFor, % Read format
     case '.grd', [Mat,Str]=Gri2Mat(str_open);
     case '.xyz', [Mat,Str]=XYZ2Mat(str_open);
     otherwise error(['Wrong data file extension: ' ReaFor '. It should  be .grd or .xyz.']);    
    end % switch ReaFor, Read format
    
    % Find blanks
    [Mat,RowGap,ColGap]=FinGap(Mat);
    [NumRow,NumCol]=size(Mat);
    % Dx: spacing in South-North direction
    % Dy: spacing in West-East direction
    Dx=(Str.Yend-Str.Yini)/(NumRow-1);
    Dy=(Str.Xend-Str.Xini)/(NumCol-1);
    
    
    % Transform DeltaT to NSS
    DeltaT2RELQ(Mat,Dx,Dy,I0_ang,D0_ang,Per,ext_para);
    
    
    FulNamR=[DirRes FilSep Nam '_R'];
    load R_mod;
    R_mod=RecMat(R_mod,FulNamR,RowGap,ColGap,Str,SavFor);
    clear R_mod;
    delete('R_mod.mat');

    FulNamE=[DirRes FilSep Nam '_E'];
    load E_mod;
    E_mod=RecMat(E_mod,FulNamE,RowGap,ColGap,Str,SavFor);
    clear E_mod;
    delete('E_mod.mat');

    FulNamL=[DirRes FilSep Nam '_L'];
    load L_mod;
    L_mod=RecMat(L_mod,FulNamL,RowGap,ColGap,Str,SavFor);
    clear L_mod;
    delete('L_mod.mat');

    FulNamQ=[DirRes FilSep Nam '_Q'];
    load Q_mod;
    Q_mod=RecMat(Q_mod,FulNamQ,RowGap,ColGap,Str,SavFor);
    clear Q_mod;
    delete('Q_mod.mat');

    % clear Haxx Haxy Haxz Hayx Hayy Hayz Zax Zay Zaz;
    
    end % function cal_mag_trans_NSS_ext
    
    
    
    
    %=======================================================================
    % FUNCTIONS 
    %=======================================================================
    
    function DeltaT2RELQ(Mat,Dx,Dy,I0_ang,D0_ang,Per,ext_para)
    [NumRow,NumCol]=size(Mat);
    % Transform DeltaT to three components
    Delta2HaxHayZa(Mat,Dx,Dy,I0_ang,D0_ang,Per,ext_para);
    % Cut the resultant matrices to original size
    load X;
    Hax=X;
    %  delete('X.mat');
    load Y;
    Hay=Y;
    %  delete('Y.mat');
    load Z;
    Za=Z;
    %  delete('Z.mat');
    %  delete('T.mat');
    delete('X.mat');
    delete('Y.mat');
    delete('Z.mat');
    clear X Y Z;
    
    % Calculate the derivates of Hax
    XGrad(Hax,Dx,Dy,Per,ext_para);
    load dataX;
    Haxx=dataX;
    
    YGrad(Hax,Dx,Dy,Per,ext_para);
    load dataY;
    Haxy=dataY;
    
    ZGrad(Hax,Dx,Dy,Per,ext_para);
    load dataZ;
    Haxz=dataZ;
    
    clear dataX dataY dataZ;
    delete('dataX.mat');
    delete('dataY.mat');
    delete('dataZ.mat');
    
    % Calculate the derivates of Hay
    Hayx=Haxy;
    
    YGrad(Hay,Dx,Dy,Per,ext_para);
    load dataY;
    Hayy=dataY;
    
    ZGrad(Hay,Dx,Dy,Per,ext_para);
    load dataZ;
    Hayz=dataZ;
    
    clear dataY dataZ;
    delete('dataY.mat');
    delete('dataZ.mat');
    
    % Calculate the derivates of Hay
    Zax=Haxz;
    
    Zay=Hayz;
    
    ZGrad(Za,Dx,Dy,Per,ext_para);
    load dataZ;
    Zaz=dataZ;
    
    clear dataZ;
    delete('dataZ.mat');
    
    % Haxx,Haxy,Haxz,Hayx,Hayy,Hayz,Zax,Zay,Zaz
    % Hax, Hay, Za
    Ta = sqrt(Hax.*Hax + Hay.*Hay + Za.*Za);
    % Calculate the R : start
    Tax = (Hax.*Haxx + Hay.*Hayx + Za.*Zax)./Ta;
    Tay = (Hax.*Haxy + Hay.*Hayy + Za.*Zay)./Ta;
    Taz = (Hax.*Haxz + Hay.*Hayz + Za.*Zaz)./Ta;
    
    R_mod = sqrt(Tax.*Tax + Tay.*Tay + Taz.*Taz);

    % Calculate the E : start
    Haxa = Haxx.*Haxx + Haxy.*Haxy + Haxz.*Haxz;
    Haya = Hayx.*Hayx + Hayy.*Hayy + Hayz.*Hayz;
    Zaa  = Zax.*Zax + Zay.*Zay + Zaz.*Zaz;
    E_mod = sqrt((Haxa + Haya + Zaa)/2);

    Q_mod = sqrt(2*E_mod.*E_mod - R_mod.*R_mod);

    L_mod = (Q_mod.*Q_mod)./Ta;


    save R_mod R_mod;
    save E_mod E_mod;
    save Q_mod Q_mod;
    save L_mod L_mod;

    % Square of X
    %T2=X.*X;
    clear R_mod E_mod Q_mod L_mod;
    clear Haxx Haxy Haxz Hayx Hayy Hayz Zax Zay Zaz;
    clear Ta Tax Tay Taz Haxa Haya Zaa;

    
    end % DeltaT2RELQ
    
    
    %=======================================================================
    function GradTensor2NSS(Haxx,Haxy,Haxz,Hayx,Hayy,Hayz,Zax,Zay,Zaz,NumRow,NumCol)
    
    for ij = 1:NumRow
        for ik = 1:NumCol
            temp_grad_tensor(1,1) = Haxx(ij,ik);
            temp_grad_tensor(1,2) = Haxy(ij,ik);
            temp_grad_tensor(1,3) = Haxz(ij,ik);
            temp_grad_tensor(2,1) = Hayx(ij,ik);
            temp_grad_tensor(2,2) = Hayy(ij,ik);
            temp_grad_tensor(2,3) = Hayz(ij,ik);
            temp_grad_tensor(3,1) = Zax(ij,ik);
            temp_grad_tensor(3,2) = Zay(ij,ik);
            temp_grad_tensor(3,3) = Zaz(ij,ik);
    
            temp_eigenvalue = eig(temp_grad_tensor);
            lamoda1 = max(temp_eigenvalue);
            lamoda3 = min(temp_eigenvalue);
            lamoda2 = sum(temp_eigenvalue) - lamoda1 - lamoda3;
    
            NSS(ij,ik) = sqrt((-1)*lamoda2*lamoda2 - lamoda1*lamoda3);
        end
    end
    save NSS NSS;
    % Square of X
    %T2=X.*X;
    clear NSS;
    disp('Recovered NSS')
    
    end % GradTensor2NSS
    
    
    %=======================================================================
    function Delta2HaxHayZa(data,Dx,Dy,Incl,Decl,Per,ext_para)
    
    [NumRow NumCol]=size(data);
    if ext_para == 1 %余弦函数扩边
        [limrow1 limrow2]=Eve10Ext(NumRow,Per);
        [limcol1 limcol2]=Eve10Ext(NumCol,Per);
        data1=ExtCos2(data,NumRow,NumCol,limrow1,limrow2,limcol1,limcol2); 
    elseif ext_para == 2 % 多项式扩边
        line = 2^(floor(log2(NumRow))+2); %扩边行数
        list = 2^(floor(log2(NumCol))+2); %扩边列数
        [data1,ext1,ext2,ext3,ext4] = polynomial_extension_2D(data,line,list,2,2,3);
    end
    %save data1 data1
    G=fft2(data1);
    [nr,nc]=size(G);
    data=[];
    data1=[];
    % Half sizes
    % nr  half
    nrh=nr/2;
    % nr  half Minus 1
    nrhM=nrh-1;
    % nr  halfPlus 1
    nrhP=nrh+1;
    nrhP2= nrh+2;
     
    % nc  half
    nch=nc/2;
    % nc  half Minus 1
    nchM=nch-1;
    % nc  half Plus 1
    nchP=nch+1;
    nchP2=nch+2;
    % Clean memory (All the upper half part)
    G(nrhP2:nr,:)=[];
     
    % Take quadrants G11 and G21 that characterize the whole spectrum G
    G11=G(:,1:nchP);
    save G11 G11;
    clear G11;
     
    G21=G(:,nchP2:nc);
    save G21 G21;
    clear G21;
     
    clear G;
    
    % Calculate the directional cosines: s,t,p
    DR=Deg2Rad(Decl);
    IR=Deg2Rad(Incl);
    s=sin(IR);
    t=cos(IR)*cos(DR);
    p=cos(IR)*sin(DR);
    
    % Calculate the elementary frequency kx0, ky0
    kx0=(2*pi)/((nr-1)*Dx);
    ky0=(2*pi)/((nc-1)*Dy);
    
    %-----------------------------------------------------------------------
    % Quadrant 11
    %-----------------------------------------------------------------------
    % Calculate filters from DeltaT to X, Y, Z
    % Here FZ11 refers to K11
    % Here FX11 refers to InvB11
    [FX11,FZ11]=CalInvB(nrhP,nchP,s,t,p,kx0,ky0);
    % Z component filter
    FZ11=FZ11.*FX11; 
     
    % Calculation of equation (21), quadrant 11
    %----------------------------------------------------------------------
    % Calculate the Fourier transform of the Z component
    load G11;
    FZ11=G11.*FZ11;
    save FZ11 FZ11; 
    clear FZ11; 
    
    % Calculation of equation (20), quadrant 11
    %----------------------------------------------------------------------
    % Y component filter
    FY11=CalFy(nchP,FX11,ky0);
    % Calculate the Fourier transform of Y component
    FY11=G11.*FY11;
    save FY11 FY11;
    clear FY11; 
    
    % Calculation of equation (19), quadrant 11
    %----------------------------------------------------------------------
    % X component filter
    FX11=CalFx(nrhP,FX11,kx0);
    % Calculate the Fourier transforms of X component
    FX11=G11.*FX11;
    save FX11 FX11;
    clear FX11; 
    clear G11;
    delete('G11.mat');
     
    disp('Finished quadrant 11')
     
    %----------------------------------------------------------------------- 
    % Quadrant 21
    %-----------------------------------------------------------------------
    % Calculate filter from DeltaT to X, Y, Z 
    % Here FZ21 refers to K21
    % Here FX21 refers to InvB21
    
    [FX21,FZ21]=CalInvBR(nrhP,nchM,s,t,p,kx0,ky0);
    % Z component filter 
    FZ21=FZ21.*FX21; 
    
    % Calculation of equation (21), quadrant 21
    %----------------------------------------------------------------------
    % Calculate the Fourier transform of Z component
    load G21;
    FZ21=G21.*FZ21;
    save FZ21 FZ21; 
    clear FZ21; 
    
    % Calculation of equation (20), quadrant 21
    %----------------------------------------------------------------------
    % Y component filter
    FY21=CalFyR(FX21,ky0);
    % Calculate the Fourier transforms of Y component
    FY21=G21.*FY21;
    save FY21 FY21;
    clear FY21; 
     
    % Calculation of equation (19), quadrant 21
    %----------------------------------------------------------------------
    % X component of the filter
    FX21=CalFxR(FX21,kx0);
    % Calculate the Fourier Transforms of the X component
    FX21=G21.*FX21;
     
    save FX21 FX21;
    clear FX21; 
    clear G21;
    delete('G21.mat');
     
    disp('Finished quadrant 21') 
    
    %-----------------------------------------------------------------------
    % Fourier transforms of Quadrant 12 Upper Left (symmetries)
    %-----------------------------------------------------------------------
    % Fourier transform of Z component
    load FZ21;
      
    FZ12=FZ21(2:nrh,:);
    clear FZ21;
    FZ12=flipud(FZ12);
    FZ12=fliplr(FZ12);
    FZ12=conj(FZ12);
     
    load FZ11;
    FZ02=FZ11(2:nrh,1);
    FZrc2=FZ11(2:nrh,nchP);
    clear FZ11;
     
    FZ02=flipud(FZ02);
    FZ02=conj(FZ02);
     
    FZrc2=flipud(FZrc2);
    FZrc2=conj(FZrc2);
      
    FZ12=[FZ02 FZ12 FZrc2];
     
    save FZ12 FZ12;
    clear FZ12;
     
    % Fourier transform of X component
    load FX21;
    FX12=FX21(2:nrh,:);
    clear FX21;
    FX12=flipud(FX12);
    FX12=fliplr(FX12);
    FX12=conj(FX12);
     
    load FX11;
    FX02=FX11(2:nrh,1);
    FXrc2=FX11(2:nrh,nchP);
    clear FX11;
     
    FX02=flipud(FX02);
    FX02=conj(FX02);
     
    FXrc2=flipud(FXrc2);
    FXrc2=conj(FXrc2);
      
    FX12=[FX02 FX12 FXrc2];
     
    save FX12 FX12;
    clear FX12;
     
    % Fourier transform of Y component
    load FY21;
    FY12=FY21(2:nrh,:);
    clear FY21;
    FY12=flipud(FY12);
    FY12=fliplr(FY12);
    FY12=conj(FY12);
     
    load FY11;
    FY02=FY11(2:nrh,1);
    FYrc2=FY11(2:nrh,nchP);
    clear FY11;
     
    FY02=flipud(FY02);
    FY02=conj(FY02);
     
    FYrc2=flipud(FYrc2);
    FYrc2=conj(FYrc2);
     
    FY12=[FY02 FY12 FYrc2];
    save FY12 FY12;
    clear FY12;
     
    disp('Finished quadrant 12')
    %-----------------------------------------------------------------------
    % Fourier transform of Z component
    load FZ11;
    FZ22=FZ11(2:nrh,2:nch);
    clear FZ11;
    FZ22=flipud(FZ22);
    FZ22=fliplr(FZ22);
    FZ22=conj(FZ22);
    save FZ22 FZ22;
    clear FZ22;
     
    % Fourier transform of X component
    load FX11;
    FX22=FX11(2:nrh,2:nch);
    clear FX11;
    FX22=flipud(FX22);
    FX22=fliplr(FX22);
    FX22=conj(FX22);
    save FX22 FX22;
    clear FX22;
     
    % Fourier transform of Y component
    load FY11;
    FY22=FY11(2:nrh,2:nch);
    clear FY11;
    FY22=flipud(FY22);
    FY22=fliplr(FY22);
    FY22=conj(FY22);
    save FY22 FY22;
    clear FY22;
    disp('Finished quadrant 22')
    %-----------------------------------------------------------------------
    % Recover the Fourier transform of X component
    %-----------------------------------------------------------------------
    load FX11;
    load FX21;
    load FX12;
    load FX22;
      
    %FX=[ FX11    FX21;
    %     FX12    FX22];
    FX=zeros(nr,nc);
    FX(1:nrhP,1:nchP)=FX11;
    FX(1:nrhP,1+nchP:end)=FX21;
    FX(1+nrhP:end,1:nchP)=FX12;
    FX(1+nrhP:end,1+nchP:end)=FX22;
     
    clear FX11;  
    clear FX12;
    delete('FX11.mat');
    delete('FX12.mat');
    clear FX21;
    clear FX22;    
    delete('FX21.mat');
    delete('FX22.mat');
    
    % Calculation of equation (22), X component in the space domain
    %----------------------------------------------------------------------- 
    X1=ifft2(FX);
    X1=real(X1);
    %X=X(1:NumRow,1:NumCol);
    for ii=1:NumRow
        for jj=1:NumCol
            if ext_para == 1
                X(ii,jj)=X1(ii+limrow1,jj+limcol1);
            elseif ext_para == 2
                X(ii,jj)=X1(ii+ext4,jj+ext1);
            end
        end
    end
    clear X1;    	
    save X X;
    % Square of X
    %T2=X.*X;
    clear X;
    disp('Recovered X')
    %-----------------------------------------------------------------------
    % Recover Fourier transform of Y component  
    %-----------------------------------------------------------------------
    load FY11;
    load FY21;
    load FY12;
    load FY22;
      
    %FY=[ FY11    FY21;
    %     FY12    FY22];
    FY=zeros(nr,nc);
    FY(1:nrhP,1:nchP)=FY11;
    FY(1:nrhP,1+nchP:end)=FY21;
    FY(1+nrhP:end,1:nchP)=FY12;
    FY(1+nrhP:end,1+nchP:end)=FY22;
      
    clear FY11;    
    clear FY12;
    delete('FY11.mat');
    delete('FY12.mat');
     
    clear FY21;
    clear FY22;  
    delete('FY21.mat');
    delete('FY22.mat');
    
    % Calculation of equation (23), Y component in the space domain
    %----------------------------------------------------------------------- 
    Y1=ifft2(FY);
    Y1=real(Y1);
    %Y=Y(1:NumRow,1:NumCol);
    for ii=1:NumRow
        for jj=1:NumCol
            if ext_para == 1
                Y(ii,jj)=Y1(ii+limrow1,jj+limcol1);
            elseif ext_para == 2
                Y(ii,jj)=Y1(ii+ext4,jj+ext1);
            end
        end
    end
    clear Y1;
    save Y Y;
    % Square of Y
    %Y=Y.*Y;
    %T2=T2+Y;
    
    clear Y;
    disp('Recovered Y')
    %-----------------------------------------------------------------------
    % Recover Fourier Transform of Z component
    %-----------------------------------------------------------------------
    load FZ11;
    load FZ21;
    load FZ12;
    load FZ22;
      
    % FZ=[ FZ11    FZ21;
    %     FZ12    FZ22];
    FZ=zeros(nr,nc);
    FZ(1:nrhP,1:nchP)=FZ11;
    FZ(1:nrhP,1+nchP:end)=FZ21;
    FZ(1+nrhP:end,1:nchP)=FZ12;
    FZ(1+nrhP:end,1+nchP:end)=FZ22;
      
    clear FZ11;    
    clear FZ12;
    delete('FZ11.mat');
    delete('FZ12.mat');
    
    clear FZ21;
    clear FZ22;   
    delete('FZ21.mat');   
    delete('FZ22.mat');
    
    % Calculation of equation (24), Z component in the space domain
    %----------------------------------------------------------------------
    Z1=ifft2(FZ);
    Z1=real(Z1);
    %Z=Z(1:NumRow,1:NumCol);
    for ii=1:NumRow
        for jj=1:NumCol
            if ext_para == 1
                Z(ii,jj)=Z1(ii+limrow1,jj+limcol1);
            elseif ext_para == 2
                Z(ii,jj)=Z1(ii+ext4,jj+ext1);
            end
        end
    end
    clear Z1;
    save Z Z;
    % Square of Z
    %Z=Z.*Z;
    %T2=T2+Z;
    
    clear Z;
    %T=sqrt(T2);
    disp('Recovered Z')
    %save T T;
    
    end % function Delta2HaxHayZa
    
    
    function XGrad(data,Dx,Dy,Per,ext_para)
    [NumRow NumCol]=size(data);
    if ext_para == 1 %余弦函数扩边
        [limrow1 limrow2]=Eve10Ext(NumRow,Per);
        NumRowExt = NumRow + limrow1 + limrow2;
        [limcol1 limcol2]=Eve10Ext(NumCol,Per);
        NumColExt = NumCol + limcol1 + limcol2;
        data1=ExtCos2(data,NumRow,NumCol,limrow1,limrow2,limcol1,limcol2); 
    elseif ext_para == 2 % 多项式扩边
        NumRowExt = 2^(floor(log2(NumRow))+2); %扩边行数
        NumColExt = 2^(floor(log2(NumCol))+2); %扩边列数
        [data1,ext1,ext2,ext3,ext4] = polynomial_extension_2D(data,NumRowExt,NumColExt,2,2,3);
    end
    %line = 2^(floor(log2(NumRow))+2); %扩边行数
    %list = 2^(floor(log2(NumCol))+2); %扩边列数
    %[dataExt,ext1,ext2,ext3,ext4] = polynomial_extension_2D(data,line,list,2,2,3);
    data = [];
    disp('Started calculating FFT of X grad component field');
    % Fourier transform of extended data
    dataExt=fft2(data1);
    Par1=NumRowExt/2+1;
    Par2=NumColExt/2+1;
    Par2p1=Par2+1;
    Par3=NumColExt/2-1;
    FX11=dataExt(1:Par1,1:Par2);
    FX21=dataExt(1:Par1, Par2p1:NumColExt);
    clear data1;
    
    kx0e=(2*pi)/((NumRowExt-1)*Dx);
    dataX1=FunXFull(FX11,FX21,Par1,Par2,Par3,NumRowExt,NumColExt,kx0e);
    %dataX=dataX(1:line,1:list);
    for ii=1:NumRow
        for jj=1:NumCol
            if ext_para == 1
                dataX(ii,jj)=dataX1(ii+limrow1,jj+limcol1);
            elseif ext_para == 2
                dataX(ii,jj)=dataX1(ii+ext4,jj+ext1);
            end
            
        end
    end
    clear dataX1;
    save dataX dataX; 
    clear dataX ;
    
    
    
    
    clear FX11;
    % delete('FX11.mat');
    clear FX21;
    % delete('FX21.mat');
    end % function XGrad
    
    function YGrad(data,Dx,Dy,Per,ext_para)
    [NumRow NumCol]=size(data);
    if ext_para == 1 %余弦函数扩边
        [limrow1 limrow2]=Eve10Ext(NumRow,Per);
        NumRowExt = NumRow + limrow1 + limrow2;
        [limcol1 limcol2]=Eve10Ext(NumCol,Per);
        NumColExt = NumCol + limcol1 + limcol2;
        data1=ExtCos2(data,NumRow,NumCol,limrow1,limrow2,limcol1,limcol2); 
    elseif ext_para == 2 % 多项式扩边
        NumRowExt = 2^(floor(log2(NumRow))+2); %扩边行数
        NumColExt = 2^(floor(log2(NumCol))+2); %扩边列数
        [data1,ext1,ext2,ext3,ext4] = polynomial_extension_2D(data,NumRowExt,NumColExt,2,2,3);
    end
    data = [];
    disp('Started calculating FFT of Y grad component field');
    % Fourier transform of extended data
    dataExt=fft2(data1);
    Par1=NumRowExt/2+1;
    Par2=NumColExt/2+1;
    Par2p1=Par2+1;
    Par3=NumColExt/2-1;
    FX11=dataExt(1:Par1,1:Par2);
    FX21=dataExt(1:Par1, Par2p1:NumColExt);
    clear data1;
    
    ky0e=(2*pi)/((NumColExt-1)*Dy);
    dataY1=FunYFull(FX11,FX21,Par1,Par2,Par3,NumRowExt,NumColExt,ky0e);
    %dataY=dataY(1:NumRow,1:NumCol);  
    for ii=1:NumRow
        for jj=1:NumCol
    %        dataY(ii,jj)=dataY1(ii+ext4,jj+ext1);
            if ext_para == 1
                dataY(ii,jj)=dataY1(ii+limrow1,jj+limcol1);
            elseif ext_para == 2
                dataY(ii,jj)=dataY1(ii+ext4,jj+ext1);
            end    
        end
    end
    clear dataY1;
    save dataY dataY; 
    clear dataY;
    clear FX11;
    % delete('FX11.mat');
    clear FX21;
    % delete('FX21.mat');
    
    end
    
    
    
    function ZGrad(data,Dx,Dy,Per,ext_para)
    [NumRow NumCol]=size(data);
    if ext_para == 1 %余弦函数扩边
        [limrow1 limrow2]=Eve10Ext(NumRow,Per);
        NumRowExt = NumRow + limrow1 + limrow2;
        [limcol1 limcol2]=Eve10Ext(NumCol,Per);
        NumColExt = NumCol + limcol1 + limcol2;
        data1=ExtCos2(data,NumRow,NumCol,limrow1,limrow2,limcol1,limcol2); 
    elseif ext_para == 2 % 多项式扩边
        NumRowExt = 2^(floor(log2(NumRow))+2); %扩边行数
        NumColExt = 2^(floor(log2(NumCol))+2); %扩边列数
        [data1,ext1,ext2,ext3,ext4] = polynomial_extension_2D(data,NumRowExt,NumColExt,2,2,3);
    end
    
    data = [];
    disp('Started calculating FFT of Z grad component field');
    % Fourier transform of extended data
    dataExt=fft2(data1);
    Par1=NumRowExt/2+1;
    Par2=NumColExt/2+1;
    Par2p1=Par2+1;
    Par3=NumColExt/2-1;
    FX11=dataExt(1:Par1,1:Par2);
    FX21=dataExt(1:Par1, Par2p1:NumColExt);
    clear data1;
    
    kx0e=(2*pi)/((NumRowExt-1)*Dx);
    ky0e=(2*pi)/((NumColExt-1)*Dy);
    dataZ1=FunZFull(FX11,FX21,Par1,Par2,Par3,NumRowExt,NumColExt,kx0e,ky0e);
    %dataZ1=dataZ(1:NumRow,1:NumCol);  
    for ii=1:NumRow
        for jj=1:NumCol
    %        dataZ(ii,jj)=dataZ1(ii+ext4,jj+ext1);
            if ext_para == 1
                dataZ(ii,jj)=dataZ1(ii+limrow1,jj+limcol1);
            elseif ext_para == 2
                dataZ(ii,jj)=dataZ1(ii+ext4,jj+ext1);
            end
        end
    end
    clear dataZ1;
    
    save dataZ dataZ; 
    clear dataZ;
    clear FX11;
    % delete('FX11.mat');
    clear FX21;
    % delete('FX21.mat');
    
    end
    
    
    function Fdz=CalFdz(Ff,nx,ny,kx0,ky0)
    %=======================================================================
    % Calculates Fy of the transfer function d/dy for quadrant F11 
    % (Left lower part)
    Fz=zeros(nx,ny);
    for x=1:nx
    for y=1:ny
       kx=(x-1)*kx0; 
       ky=(y-1)*ky0;   
       Fz(x,y)=sqrt(kx*kx + ky*ky);
    end % for y=1:ny
    end % for x=1:nx
    Fdz=Ff.*Fz;
    end
    
    
    function FdzR=CalFdzR(Ff,nx,ny,kx0,ky0)
    %=======================================================================
    % Calculates Fx of the transfer function d/dx for quadrant F12 
    % (Right lower part)
    Fz=zeros(nx,ny);
    for x=1:nx
    for y=1:ny  
       kx=(x-1)*kx0; 
       ky=(y-ny-1)*ky0;
       Fz(x,y)=sqrt(kx*kx + ky*ky);
    end % for y=1:ny
    end
    FdzR=Ff.*Fz;
    end
    %=======================================================================
    function Fdx=CalFdx(Ff,nx,ny,kx0)
    %=======================================================================
    % Calculates Fx of the transfer function d/dx for quadrant F11 
    % (Left lower part)
    Fx=zeros(nx,ny);
    for x=1:nx
       kx=(x-1)*kx0;   
       Fx(x,:)=i*kx;
    end %for x=1:nx
    Fdx=Ff.*Fx;
    end
    %=======================================================================
    function Fdy=CalFdy(Ff,nx,ny,ky0)
    %=======================================================================
    % Calculates Fy of the transfer function d/dy for quadrant F11 
    % (Left lower part)
    Fy=zeros(nx,ny);
    for y=1:ny
       ky=(y-1)*ky0;   
       Fy(:,y)=i*ky;
    end % for y=1:ny
    Fdy=Ff.*Fy;
    end
    %=======================================================================
    function FdyR=CalFdyR(Ff,nx,ny,ky0)
    %=======================================================================
    % Calculates Fx of the transfer function d/dx for quadrant F12 
    % (Right lower part)
    Fy=zeros(nx,ny);
    for y=1:ny  
       ky=(y-ny-1)*ky0; 
       Fy(:,y)=i*ky;
    end % for y=1:ny
    FdyR=Ff.*Fy;
    end
    %=======================================================================
    function FunX=FunXFull(FfL,FfR,nxLR,nyL,nyR,nx,ny,kx0)
    %=======================================================================
    % Calculates x horizontal derivative of a function in Fourier domain  
    % in the 4 quadrants
    nrh=nx/2;
    nch=ny/2;
    Xx11=CalFdx(FfL,nxLR,nyL,kx0);
    Xx21=CalFdx(FfR,nxLR,nyR,kx0);
    
    % Calculate upper quadrants and space domain Xx and Xy
    Xx12=Xx21(2:nrh,:);
    Xx12=flipud(Xx12);
    Xx12=fliplr(Xx12);
    Xx12=conj(Xx12);
    Xx02=Xx11(2:nrh,1);
    Xxrc2=Xx11(2:nrh,nyL);
    Xx02=flipud(Xx02);
    Xx02=conj(Xx02);
    Xxrc2=flipud(Xxrc2);
    Xxrc2=conj(Xxrc2);
    Xx12=[Xx02 Xx12 Xxrc2];
    Xx22=Xx11(2:nrh,2:nch);
    Xx22=flipud(Xx22);
    Xx22=fliplr(Xx22);
    Xx22=conj(Xx22);
    FunX=zeros(nx,ny);
    FunX(1:nxLR,1:nyL)=Xx11;
    FunX(1:nxLR,1+nyL:end)=Xx21;
    FunX(1+nxLR:end,1:nyL)=Xx12;
    FunX(1+nxLR:end,1+nyL:end)=Xx22;
    clear Xx11;    
    clear Xx12;
    clear Xx21;
    clear Xx22;   
    
    % Spectrum back to space domain
    FunX=ifft2(FunX);
    FunX=real(FunX);
    end
    %=======================================================================
    function FunY=FunYFull(FfL,FfR,nxLR,nyL,nyR,nx,ny,ky0)
    %=======================================================================
    % Calculates horizontal derivative y of a function in Fourier domain
    % in the 4 quadrants
    nrh=nx/2;
    nch=ny/2;
    Xx11=CalFdy(FfL,nxLR,nyL,ky0);
    Xx21=CalFdyR(FfR,nxLR,nyR,ky0);
    
    % Calculate upper quadrants and space domain Xx and Xy
    Xx12=Xx21(2:nrh,:);
    Xx12=flipud(Xx12);
    Xx12=fliplr(Xx12);
    Xx12=conj(Xx12);
    Xx02=Xx11(2:nrh,1);
    Xxrc2=Xx11(2:nrh,nyL);
    Xx02=flipud(Xx02);
    Xx02=conj(Xx02);
    Xxrc2=flipud(Xxrc2);
    Xxrc2=conj(Xxrc2);
    Xx12=[Xx02 Xx12 Xxrc2];
    Xx22=Xx11(2:nrh,2:nch);
    Xx22=flipud(Xx22);
    Xx22=fliplr(Xx22);
    Xx22=conj(Xx22);
    FunY=zeros(nx,ny);
    FunY(1:nxLR,1:nyL)=Xx11;
    FunY(1:nxLR,1+nyL:end)=Xx21;
    FunY(1+nxLR:end,1:nyL)=Xx12;
    FunY(1+nxLR:end,1+nyL:end)=Xx22;
    clear Xx11;    
    clear Xx12;
    clear Xx21;
    clear Xx22;  
    
    % Xx in the space domain
    FunY=ifft2(FunY);
    FunY=real(FunY);
    end
    
    
    function FunZ=FunZFull(FfL,FfR,nxLR,nyL,nyR,nx,ny,kx0,ky0)
    %=======================================================================
    % Calculates x horizontal derivative of a function in Fourier domain  
    % in the 4 quadrants
    nrh=nx/2;
    nch=ny/2;
    Xx11=CalFdz(FfL,nxLR,nyL,kx0,ky0);
    Xx21=CalFdzR(FfR,nxLR,nyR,kx0,ky0);
    
    % Calculate upper quadrants and space domain Xx and Xy
    Xx12=Xx21(2:nrh,:);
    Xx12=flipud(Xx12);
    Xx12=fliplr(Xx12);
    Xx12=conj(Xx12);
    Xx02=Xx11(2:nrh,1);
    Xxrc2=Xx11(2:nrh,nyL);
    Xx02=flipud(Xx02);
    Xx02=conj(Xx02);
    Xxrc2=flipud(Xxrc2);
    Xxrc2=conj(Xxrc2);
    Xx12=[Xx02 Xx12 Xxrc2];
    Xx22=Xx11(2:nrh,2:nch);
    Xx22=flipud(Xx22);
    Xx22=fliplr(Xx22);
    Xx22=conj(Xx22);
    FunZ=zeros(nx,ny);
    FunZ(1:nxLR,1:nyL)=Xx11;
    FunZ(1:nxLR,1+nyL:end)=Xx21;
    FunZ(1+nxLR:end,1:nyL)=Xx12;
    FunZ(1+nxLR:end,1+nyL:end)=Xx22;
    clear Xx11;    
    clear Xx12;
    clear Xx21;
    clear Xx22;   
    
    % Spectrum back to space domain
    FunZ=ifft2(FunZ);
    FunZ=real(FunZ);
    end
    %======================================================================
    
    
    function Mat=RecMat(Map,NamFil,I,J,Str,SavFor)
    %=======================================================================
    % Recover original size
    % and postprocessing the mat files
    
    Map=Map(1:Str.NumRow,1:Str.NumCol);
    
    % Calculate the extrema before the recovering the blanks
    Min=min(min(Map));
    Max=max(max(Map));
    
    Str.Min=Min;
    Str.Max=Max;
    
    % Recover the blanks
    % Fills the gaps with 1.70141E+038  
    ValGap=1.70141E+038;
    Mat=RecGap(Map,ValGap,I,J);
    
    % Save matrix data
    NamFilExt=[NamFil '.' SavFor];
    switch SavFor,
     case 'grd', % GS ASCII grid format     
      Mat2Gri(NamFilExt,Mat,Str);
     case 'xyz', % XYZ format  
      Mat2XYZ(NamFilExt,Mat,Str);
     case '',    % Default 
     otherwise error(['Wrong output file extension: ' SavFor '. It should be grd or xyz.']);    
    end %switch SavFor,
    
    % Postprocessing the mat files
    % Fill the gaps of the mat files with NaN 
    ValGap=NaN;
    Mat=RecGap(Map,ValGap,I,J);
    end
    %=======================================================================
    
    
    function gExt=ExtCos(g,NumRow,NumCol,NumRowExt,NumColExt)
    %=======================================================================
    % Extension of the grid with half a cosine function
    % Cosine taper right and upper side
    s=0.5*(g(:,NumCol)+g(:,1));
    %r=0.5*(g(:,NumCol)-g(:,1));
    lim=(NumColExt-NumCol)/2;
    rlamb=lim+1;
    arg=pi/rlamb;
    for k=1:lim
     ge=s+r*cos(arg*k);
     g(:,NumCol+k)=ge;
    end %for k=1:lim
    g=g';
    % Cosine taper upper side
    s=0.5*(g(:,NumRow)+g(:,1));
    r=0.5*(g(:,NumRow)-g(:,1));
    lim=NumRowExt-NumRow;
    rlamb=lim+1;
    arg=pi/rlamb;
    for k=1:lim
     ge=s+r*cos(arg*k);
     g(:,NumRow+k)=ge;
    end %for k=1:lim
    gExt=g';
    end % function gExt=ExtCos
    
    function gExt=ExtCos2(g,NumRow,NumCol,limrow1,limrow2,limcol1,limcol2)
    %=======================================================================
    % Extension of the grid with half a cosine function
    % Cosine taper right and upper side
    
    %r=0.5*(g(:,NumCol)-g(:,1));
    
    rlamb=limcol1+1;
    arg=pi/rlamb;
    
    sl=0.5*(g(:,NumCol)+g(:,1));
    s=0.5*(g(:,1)+sl(:,1));
    r=0.5*(g(:,1)-sl(:,1));
    for k = 1:limcol1
        ge = s-r*cos(arg*k);
        g_temp(:,k) = ge;
    end
    g_temp(:,limcol1+1:limcol1+NumCol) = g;
    
    rlamb=limcol2+1;
    arg=pi/rlamb;
    s=0.5*(g(:,NumCol)+sl(:,1));
    r=0.5*(g(:,NumCol)-sl(:,1));
    for k = 1:limcol2
        ge = s+r*cos(arg*k);
        g_temp(:,NumCol+limcol1+k) = ge;
    end
    
    g_temp1 = g_temp';
    clear g_temp ge;
    
    
    
    
    rlamb=limrow1+1;
    arg=pi/rlamb;
    sl=0.5*(g_temp1(:,NumRow)+g_temp1(:,1));
    s=0.5*(g_temp1(:,1)+sl(:,1));
    r=0.5*(g_temp1(:,1)-sl(:,1));
    for k = 1:limrow1
        ge = s-r*cos(arg*k);
        g_temp(:,k) = ge;
    end
    g_temp(:,limrow1+1:limrow1+NumRow) = g_temp1;
    
    rlamb=limrow2+1;
    arg=pi/rlamb;
    s=0.5*(g_temp1(:,NumRow)+sl(:,1));
    r=0.5*(g_temp1(:,NumRow)-sl(:,1));
    for k = 1:limrow2
        ge = s+r*cos(arg*k);
        g_temp(:,NumRow+limrow1+k) = ge;
    end
    
    gExt=g_temp';
    end % function gExt=ExtCos2
    
    
    %=======================================================================
    
    %=======================================================================
    function [Mat,I,J]=FinGap(Mat)
    %=======================================================================
    % Find gaps in a matrix (from GS ASCII grid format or XYZ format)  
    
    % OUTPUTS:
    % MatFil : Matrix with filled blanks
    % MatGap : Matrix with the positions of the blanks
    
    ValGap=1.70141E+038;
    
    FilGap=0;
    
    % Find blanks
    [I,J]=find(Mat==ValGap);
    
    % Fill the blanks
    % Mat(I,J)=FilGap;
    NumGap=size(I,1);
    for k=1:NumGap
     i=I(k);
     j=J(k);
     Mat(i,j)=FilGap;
    end % for i=1:NumGap
    end % function [Mat,I,J]=FinGap(Mat)
    %======================================================================
    function [Mat,Str]=Gri2Mat(NamFil)
    %======================================================================
    % Covert GS ASCII grid to Matlab file format
    % OUTPUT
    % Mat: Matrix with the data
    % Str: Structure of data
    
    % Open the file
    pF=fopen(NamFil,'r');
    
    % Read DSAA label 
    Lin=fgetl(pF);
    Tok=strtok(Lin);
    if ~strcmp('DSAA',Tok)
     error(['The data file, ' NamFil ' is not in ASCII Grid File Format']);   
    end %if ~strcmp('DSAA',Tok)
    
    % Read number of columns and rows
    Lin=fgetl(pF);
    [StrCol,StrRow]=strtok(Lin);
    NumRow=str2num(StrRow);
    NumCol=str2num(StrCol);
    
    % Read the X coordinates (West-East)
    Lin=fgetl(pF);
    [StrXini,StrXend]=strtok(Lin);
    Xini=str2num(StrXini);
    Xend=str2num(StrXend);
    
    % Read the Y coordinates (South-North)
    Lin=fgetl(pF);
    [StrYini,StrYend]=strtok(Lin);
    Yini=str2num(StrYini);
    Yend=str2num(StrYend);
    
    % Read the Field Extrema
    Lin=fgetl(pF);
    [StrMin,StrMax]=strtok(Lin);
    Min=str2num(StrMin);
    Max=str2num(StrMax);
    
    % Read the data
    Mat=fscanf(pF,'%f',[NumCol,NumRow]);
    
    % Transpose 
    Mat=Mat';
    
    % Pack structure
    Str.NumRow=NumRow;
    Str.NumCol=NumCol;
    
    Str.Xini=Xini;
    Str.Xend=Xend;
    
    Str.Yini=Yini;
    Str.Yend=Yend;
    
    Str.Min=Min;
    Str.Max=Max;
    
    % Close file -----------------------------------------------------------
    fclose(pF);
    end % function [Mat,Str]=Gri2Mat(NamFil)
    %=======================================================================
    function Mat2Gri(NamFil,Mat,Str)
    %=======================================================================
    % Convert Matlab to GS ASCII grid file format
    
    % INPUT
    % NamFil: Name of the grid file 
    % Mat:    Matrix with the data
    % Str:    Structure of data
    
    % DSAA format 
    % ASCII grid files [.GRD] contain five header lines that provide 
    % information about the size and limits of the grid, followed by a list
    % of Z values. 
    % The fields within ASCII grid files must be space delimited.
    %
    % The listing of Z values follows the header information in the file. 
    % The Z values are stored in row-major order starting with the minimum 
    % Y coordinate. 
    % The first Z value in the grid file corresponds to the lower left corner
    % of the map. 
    % This can also be thought of as the southwest corner of the map, or, 
    % more specifically, the grid node of minimum X and minimum Y. 
    % The second Z value is the next adjacent grid node in the same row 
    % (the same Y coordinate but the next higher X coordinate).
    % When the maximum X value is reached in the row,
    % the list of Z values continues with the next higher row, 
    % until all the rows of Z values have been included. 
    %
    % The general format of an ASCII grid file is:
    %
    % id The identification string DSAA that identifies the file as an ASCII
    % grid file.
    % nx ny nx is the integer number of grid lines along the X axis (columns)
    % ny is the integer number of grid lines along the Y axis (rows)
    % xlo xhi xlo is the minimum X value of the grid
    % xhi is the maximum X value of the grid
    % ylo yhi ylo is the minimum Y value of the grid
    % yhi is the maximum Y value of the grid
    % zlo zhi zlo is the minimum Z value of the grid
    % zhi is the maximum Z value of the grid
    % grid row 1
    % grid row 2
    % gridrow 3
    % these are the rows of Z values of the grid, organized in row order. 
    % Each row has a constant Y coordinate.  
    % Grid row 1 corresponds to ylo and the last grid row corresponds to yhi. 
    % Within each row, the Z values are arranged from xlo to xhi.
    
    % WARNING
    % Every physical line is writen along several lines in the file of a 
    % maximum of 10 columns
    % The next line is separated by an empty line in the file
    
    % Maximum number of columns per row
    MaxCol=10;
    
    % Unpack structure
    Xini=Str.Xini;
    Xend=Str.Xend;
    
    Yini=Str.Yini;
    Yend=Str.Yend;
    
    Min=Str.Min;
    Max=Str.Max;
    
    NumRow=Str.NumRow;
    NumCol=Str.NumCol;
    
    % Open the file
    pF=fopen(NamFil,'w');
    
    % Write DSAA label
    fprintf(pF,'DSAA\n');
    
    % Write number of columns and rows
    fprintf(pF,'%i %i\n',NumCol,NumRow);
    
    % Write the X coordinates (West-East)
    fprintf(pF,'%G %G\n',Xini,Xend);
    
    % Write the Y coordinates (South-North)
    fprintf(pF,'%G %G\n',Yini,Yend);
    
    % Write the Field Extrema
    fprintf(pF,'%G %G\n',Min,Max);
    
    % Prepare the parameters for writing the data
    
    % Calculate the number of file rows (nr) of every physical row
    if NumCol<=MaxCol
     nr=1;
    else % NumCol>MaxCol
     nr=ceil(NumCol/MaxCol);   
     nclr=mod(NumCol,MaxCol); % Number of columns of the last row  
     
     if nclr==0
      % If the number of columns is a multiple of MaxCol
      % we increment artificially the number or rows in every physical row,
      % since the number of colums in the las row is cero (is missing) 
      % in this case 
      nr=nr+1;
     end % if nclr==0   
    end % if NumCol<=10
    
    % Write data
    for r=1:NumRow
     c=0;   
     for j=1:nr-1
      % Full columns
      for k=1:MaxCol
       c=c+1;   
       fprintf(pF,'%G ',Mat(r,c));
      end % for k=1:MaxCol
      fprintf(pF,'\n');
     end % for j=1:nr-1  
     % Last row, with not full columns
     for k=1:nclr
      c=c+1;   
      fprintf(pF,'%G ',Mat(r,c));
     end %for k=1:MaxCol
     if nclr~=0
      fprintf(pF,'\n');
     end %if nclr~=0 
     % New line separating every group of lines
     fprintf(pF,'\n');
    end %for i=1:NumRow    
    
    % Close file -----------------------------------------------------------
    fclose(pF);
    end % function Mat2Gri(NamFil,Mat,Str)
    
    %=======================================================================
    function Mat2XYZ(NamFil,Mat,Str)
    %=======================================================================
    % Convert Matlab to XYZ file format (without blanks)
    
    % INPUT
    % NamFil: Name of the grid file 
    % Mat:    Matrix with the data
    % Str:    Structure of data
    
    ValGap=1.70141E+038;
    
    % Unpack structure
    Xini=Str.Xini;
    Xend=Str.Xend;
    
    Yini=Str.Yini;
    Yend=Str.Yend;
    
    Min=Str.Min;
    Max=Str.Max;
    
    NumRow=Str.NumRow;
    NumCol=Str.NumCol;
    
    % Open the file
    pF=fopen(NamFil,'w');
    
    % Increments
    DelX=(Xend-Xini)/(NumCol-1);
    DelY=(Yend-Yini)/(NumRow-1);
    
    % Write data
    for i=1:NumRow
     y=Yini+(i-1)*DelY;   
     for j=1:NumCol
      z=Mat(i,j);
      if z~=ValGap
       x=Xini+(j-1)*DelX;
       fprintf(pF,'%G %G %G\n',x,y,z);
      end % if z~=ValGap 
     end % for j=1:NumCol
    end % for i=1:NumRow    
    
    % Close file -----------------------------------------------------------
    fclose(pF);
    end 
    %=======================================================================
    function [Mat]=RecGap(Mat,ValGap,I,J)
    %=======================================================================
    % Recovers gaps with the blank values given in the variable ValGap
    
    % Fill the gaps
    NumGap=size(I,1);
    for k=1:NumGap
     i=I(k);
     j=J(k);
     Mat(i,j)=ValGap;
    end % for i=1:NumGap
    end
    %=======================================================================
    function [Mat,Str]=XYZ2Mat(NamFil)
    %=======================================================================
    % Converts XYZ to Matlab file format
    
    % OUTPUT
    % Mat: Matrix with the data
    % Str: Structure of data
    
    % Open  the file
    pF=fopen(NamFil,'r');
    
    % Read the files
    [Mat3,Count]=fscanf(pF,'%G',[3 inf]);
    % Close file -----------------------------------------------------------
    fclose(pF);
    
    Mat3=Mat3';
    
    [NumRow0 NumCol0]=size(Mat3);
           
    if NumCol0~=3
     MesErr=['The data file: ' NamFil ' has ' num2str(NumCol0) ', instead of the required 3 columns for the XYZ format'];   
     error(MesErr);   
    end % if NumCol0~=3
    
    % Look for the number of columns 
    Fin='n';
    
    Yini=Mat3(1,2);
    i=0;
    NumCol=[];
    while Fin=='n'
     i=i+1;
     yi=Mat3(i,2);
     if yi~=Yini
      Fin='y';
      NumCol=i-1; % minus 1 because I discover in the next one   
     else % yi==Yini
      if i==NumRow0
       Fin='y';
       % Strange case, vector matrix 
       NumCol=NumRow0;
      end % if i==NumRow0    
     end % if y1~=yi   
    end % while Fin=='n'
    
    % Calculate the number of rows
    NumRow=NumRow0/NumCol;
    if NumRow>floor(NumRow)
      MesErr=['The data file: ' NamFil ', does not contain the regular grid required for the XYZ format'];  
      error(MesErr);   
    end % if NumRow>ceil(NumRow)   
    
    % Calculate the rest of parameters
    Xini=Mat3(1,1);
    Yini=Mat3(1,2);
    
    Xend=Mat3(end,1);
    Yend=Mat3(end,2);
    
    Min=min(Mat3(:,3));
    Max=max(Mat3(:,3));
    
    % Build the matrix
    k=0;
    for i=1:NumRow
     for j=1:NumCol
      k=k+1;
      Mat(i,j)=Mat3(k,3);
     end % for j=1:NumRow    
    end % for i=1:NumRow    
    
    % Pack structure
    Str.NumRow=NumRow;
    Str.NumCol=NumCol;
    
    Str.Xini=Xini;
    Str.Xend=Xend;
    
    Str.Yini=Yini;
    Str.Yend=Yend;
    
    Str.Min=Min;
    Str.Max=Max;
    end
    
    
    
    %=======================================================================
    function [InvB,K]=CalInvB(nx,ny,s,t,p,kx0,ky0)
    %=======================================================================
    % Calculates InvB and K for quadrant F11 (Left lower part)
    K=zeros(nx,ny);
    B=zeros(nx,ny);
    for x=1:nx
        kx=(x-1)*kx0;
        kx2=kx*kx;
       for y=1:ny
          ky=(y-1)*ky0;
          k2=kx2+ky*ky;
          k=sqrt(k2); 
          K(x,y)=k;
          B(x,y)=k*s+i*(kx*t+ky*p);
       end % for y=1:ny    
    end % for x=1:nx   
    B(1,1)=1;
    InvB=1./B; 
    InvB(1,1)=0+eps; % W(0,0)=0!
    end
    %=======================================================================
    function Fx=CalFx(nx,InvB,kx0)
    %=======================================================================
    % Calculates Fx for quadrant F11 (Left lower part)
    ny=size(InvB,2);
    Fx=zeros(nx,ny);
    for x=1:nx
       kx=(x-1)*kx0;   
       Fx(x,:)=(i*kx)*InvB(x,:);
    end % for x=1:nx
    end
    %=======================================================================
    function Fy=CalFy(ny,InvB,ky0)
    %=======================================================================
    % Calculates Fy for quadrant F11 (Left lower part)
    nx=size(InvB,1);
    Fy=zeros(nx,ny);
    for y=1:ny
       ky=(y-1)*ky0;   
       Fy(:,y)=(i*ky)*InvB(:,y);
    end % for y=1:ny
    end
    %=======================================================================
    function [InvB,K]=CalInvBR(nx,ny,s,t,p,kx0,ky0)
    %=======================================================================
    % Calculates InvB and K for quadrant F21 (Right lower part)
    K=zeros(nx,ny);
    B=zeros(nx,ny);
    for x=1:nx
       kx=(x-1)*kx0;
       kx2=kx*kx;
       for y=1:ny
          %ky=-(ny+1-y)*ky0;
          ky=(y-ny-1)*ky0;
          k2=kx2+ky*ky;
          k=sqrt(k2); 
          K(x,y)=k;
          B(x,y)=k*s+i*(kx*t+ky*p);
        end % for y=1:ny    
    end % for x=1:nx   
    B(1,1)=1;
    InvB=1./B; 
    InvB(1,1)=0+eps;
    end
    %=======================================================================
    function Fy=CalFyR(InvB,ky0)
    %=======================================================================
    % Calculates Fy for quadrant F21 (Right lower part)
    [nx,ny]=size(InvB);
    Fy=zeros(nx,ny);
    for y=1:ny  
       ky=(y-ny-1)*ky0; 
       Fy(:,y)=(i*ky)*InvB(:,y);
    end % for y=1:ny
    end
    %=======================================================================
    function Fx=CalFxR(InvB,kx0)
    %=======================================================================
    % Calculates Fx for quadrant F21 (Right lower part)
    [nx,ny]=size(InvB);
    Fx=zeros(nx,ny);
    for x=1:nx
       kx=(x-1)*kx0;   
       Fx(x,:)=(i*kx)*InvB(x,:);
    end % for x=1:nx
    end
    %=======================================================================
    function Rad=Deg2Rad(Deg)
    %=======================================================================
    % Converts the input angle from degrees (Deg) to radians (Rad)
    Rad=pi*Deg/180;
    end
    function [ext1 ext2]=Eve10Ext(Num,Per)
    %======================================================================
    % Even 10% extension of the number Num
    Per=Per/100;
    NumExt=round(Per*Num);
    if rem(Num,2) == 1
        ext1 = NumExt;
        ext2 = NumExt+1;
    else
        ext1 = NumExt;
        ext2 = NumExt;
    end
    
    end
    
    
    function [T,k1,k2,k3,k4]=polynomial_extension_2D(t,H,L,op2,op,n)
    %函数功能：局部多项式扩边，一维或二维 
    %一维扩边顺序：先左后右；二维扩边顺序：先左后右，再上后下
    %函数返回值：扩边后数据（矩阵形式），k1、k2、k3、k4为 左右 上下的扩展的点数
    %函数参数：t待扩边数据；H、L扩边后行列数，一维是只需要H；op2=1计算一维
    %op2=2计算二维；n多项式阶数
    %op=1端点值为0，op=2端点值为原数据端点值平均值
    %[T,k1,k2,k3,k4]=kuobian_duoxiangshi(t,H,L,op2,op,n)
    
    if op2==1
        m=numel(t);
        M=H;
        T=zeros(M,1);
    
        k1=floor((M-m)/2);
        k2=ceil((M-m)/2);
    
        K1=ones(n+1,n+1);
        for i=2:n+1
            for j=2:n+1
                K1(i,j)=(k1+i-1)^(j-1);
            end
        end
    
        G1=zeros(n+1,1);
        if op==1
            G1(1)=0;
            else if op==2
                 G1(1)=(t(1)+t(m))/2;
                end
        end
    
        for i=2:n+1
            G1(i)=t(i-1);
        end
        A1=K1\G1;
    
        T1=ones(k1,n+1);
        for i=2:k1
            for j=2:n+1
                T1(i,j)=i^(j-1);
            end
        end
        T2=T1*A1;
        for i=1:k1
            T(i)=T2(i);
        end
    
        K2=ones(n+1,n+1);
        for i=2:n+1
            for j=2:n+1
                K2(i,j)=(k2+i-1)^(j-1);
            end
        end
        
        G2=zeros(n+1,1);
        if op==1
            G2(1)=0;
        else if op==2
                G2(1)=(t(1)+t(m))/2;
            end
        end
    
        for i=2:n+1
            G2(i)=t(m-i+2);
        end
        A2=K2\G2;
    
        T3=ones(k2,n+1);
        for i=2:k2
            for j=2:n+1
                T3(i,j)=i^(j-1);
            end
        end
        T4=T3*A2;
        for i=1:k2
            T(i+k1+m)=T4(k2-i+1);
        end
    
        for i=k1+1:k1+m
            T(i,1)=t(i-k1);
        end
    else if op2==2
            sizet=size(t);
    ht=sizet(1,1);
    lt=sizet(1,2);
    Xn=L;%扩边到2^n个数据
    Yn=H;
    
    T=zeros(Yn,Xn);
    
    k1=floor((Xn-lt)/2);
    k2=ceil((Xn-lt)/2);
    k3=floor((Yn-ht)/2);
    k4=ceil((Yn-ht)/2);
    
    for h=1:ht %向左扩边
        K1=ones(n+1,n+1);
        for i=2:n+1
            for j=2:n+1
                K1(i,j)=(k1+i-1)^(j-1);
            end
        end
        
        G1=zeros(n+1,1);
        if op==1
            G1(1)=0;
        else if op==2
                G1(1)=(t(h,1)+t(h,lt))/2;
            end
        end
        
        for i=2:n+1
            G1(i)=t(h,i-1);
        end
        A1=K1\G1;
        T1=ones(k1,n+1);
        for i=2:k1
            for j=2:n+1
                T1(i,j)=i^(j-1);
            end
        end
        T2=T1*A1;
        for i=1:k1
            T(h+k4,i)=T2(i);
        end
    end
    
    for h=1:ht %向右扩边
        K2=ones(n+1,n+1);
        for i=2:n+1
            for j=2:n+1
                K2(i,j)=(k2+i-1)^(j-1);
            end
        end
        
        G2=zeros(n+1,1);
        if op==1
            G2(1)=0;
        else if op==2
                G2(1)=(t(h,1)+t(h,lt))/2;
            end
        end
        
        for i=2:n+1
            G2(i)=t(h,lt-i+2);
        end
        A2=K2\G2;
        T3=ones(k2,n+1);
        for i=2:k2
            for j=2:n+1
                T3(i,j)=i^(j-1);
            end
        end
        T4=T3*A2;
        for i=1:k2
            T(h+k4,i+k1+lt)=T4(k2-i+1);
        end
    end
    
    for h=1:ht
        for l=1:lt
            T(h+k4,l+k1)=t(h,l);
        end
    end
    
    for l=1:Xn %向上扩边
        K3=ones(n+1,n+1);
        for i=2:n+1
            for j=2:n+1
                K3(i,j)=(k4+i-1)^(j-1);
            end
        end
        
        G3=zeros(n+1,1);
        if op==0
            G3(1)=0;
        else if op==2
                G3(1)=(T(1+k4,l)+T(ht+k4,l))/2;
            end
        end
        
        for i=2:n+1
            G3(i)=T(k4+i-1,l);
        end
        A3=K3\G3;
        T5=ones(k4,n+1);
        for i=2:k4
            for j=2:n+1
                T5(i,j)=i^(j-1);
            end
        end
        T6=T5*A3;
        for i=1:k4
            T(i,l)=T6(i);
        end
    end
    
    for l=1:Xn %向下扩边
        K4=ones(n+1,n+1);
        for i=2:n+1
            for j=2:n+1
                K4(i,j)=(k3+i-1)^(j-1);
            end
        end
        
        G4=zeros(n+1,1);
        if op==0
            G4(1)=0;
        else if op==2
                G4(1)=(T(1+k4,l)+T(ht+k4,l))/2;
            end
        end
        
        for i=2:n+1
            G4(i)=T(k4+ht-i+2,l);
        end
        A4=K4\G4;
        T7=ones(k3,n+1);
        for i=2:k3
            for j=2:n+1
                T7(i,j)=i^(j-1);
            end
        end
        T8=T7*A4;
        for i=1:k3
            T(i+k4+ht,l)=T8(k3-i+1);
        end
    end
        end
    end
    
    end % function polynomial_extension_2D