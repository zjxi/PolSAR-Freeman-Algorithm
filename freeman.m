%%freeman分解
clc;
clear all;
close all;
row =900;     col =1024;     data = zeros(row,col,9);

%%读入相干矩阵[T]
fIn = fopen( 'T3\T11.bin','r');
data(:,:,1) = fread(fIn,[col,row],'float').';     fclose(fIn);
fIn = fopen( 'T3\T22.bin','r');
data(:,:,2) = fread(fIn,[col,row],'float').';     fclose(fIn);
fIn = fopen('T3\T33.bin','r');
data(:,:,3) = fread(fIn,[col,row],'float').';     fclose(fIn);
fIn = fopen('T3\T12_real.bin','r');
data(:,:,4) = fread(fIn,[col,row],'float').';   fclose(fIn);
fIn = fopen( 'T3\T13_real.bin','r');
data(:,:,5) = fread(fIn,[col,row],'float').';   fclose(fIn);
fIn = fopen( 'T3\T23_real.bin','r');
data(:,:,6) = fread(fIn,[col,row],'float').';   fclose(fIn);
fIn = fopen('T3\T12_imag.bin','r');
data(:,:,7) = fread(fIn,[col,row],'float').';   fclose(fIn);
fIn = fopen('T3\T13_imag.bin','r');
data(:,:,8) = fread(fIn,[col,row],'float').';   fclose(fIn);
fIn = fopen( 'T3\T23_imag.bin','r');
data(:,:,9) = fread(fIn,[col,row],'float').';   fclose(fIn); 
clear fIn;

%%得到相干矩阵[T]
t11=data(:,:,1);
t22=data(:,:,2);
t33=data(:,:,3);
t12=data(:,:,4)+1i*data(:,:,7);
t13=data(:,:,5)+1i*data(:,:,8);
t23=data(:,:,6)+1i*data(:,:,9);   

A=(1/sqrt(2))*[1 0 1;1 0 -1;0 sqrt(2) 0];    %相干矩阵[T]和协方差矩阵的转化系数
for ii=1:900
    for jj=1:1024
        t=[t11(ii,jj) t12(ii,jj) t13(ii,jj);t12(ii,jj)' t22(ii,jj) t23(ii,jj);t13(ii,jj)' t23(ii,jj)' t33(ii,jj)];
        c=inv(A)*t*inv(A');      %%得到协方差矩阵
        fv=(3/2)*c(2,2);
        pv=(8/3)*fv;     %求pv
        if real(c(1,3))>=0
            alph=-1;
            beta=1+(c(1,1)-c(3,3))/(c(1,3)+c(3,3)-2*c(2,2));
% be=(c(1,1)+c(1,3)-2*c(2,2))/(c(1,3)+c(3,3)-2*c(2,2));


%%求ps和pd
           fs = (c(1,3) + c(3,3) - 2*c(2,2))/(1 + beta);
            ps=fs*(1+(abs(beta))^2);
            fd=c(3,3)-fs-fv;
            pd=2*fd;
        else
            beta=1;
            alph=(c(3,3)-c(1,1))/(c(3,3)-c(1,3)-c(2,2))-1;
            fd=(c(3,3)-c(1,3)-c(2,2))/(1-alph);
            fs=c(3,3)-fd-fv;
            pd=fd*(1+(abs(alph))^2);
            ps=2*fs;
        end
%%得到Pv、Ps、Pd三个分量


        Pv(ii,jj)=pv;
        Ps(ii,jj)=ps;
        Pd(ii,jj)=pd;  
    end
end
freem(:,:,1)=Pd;
freem(:,:,2)=Pv;
freem(:,:,3)=Ps;
figure;imshow(Pv);title('Pv');
figure;imshow(Pd);title('Pd');
figure;imshow(Ps);title('Ps');
figure;imshow(freem);title('freeman-decomposition');  %混色后的RGB图像
