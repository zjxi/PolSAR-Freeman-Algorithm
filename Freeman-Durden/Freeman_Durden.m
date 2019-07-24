% Freeman-Durden分类
clc;
clear all;
close all;
height=900;
width=1024;

% matlab文件操作 
file_id=fopen('T3\T11.bin','r');   
status=fseek(file_id,0,'bof');   %定位文件位置指针
for ii=1:height
    pix=fread(file_id,width,'float32');
    T11(ii,:)=pix;
end
fclose(file_id);

file_id=fopen('T3\T22.bin','r');
status=fseek(file_id,0,'bof');
for ii=1:height
    pix=fread(file_id,width,'float32');
    T22(ii,:)=pix;
end
fclose(file_id);

file_id=fopen('T3\T13_real.bin','r');
status=fseek(file_id,0,'bof');
for ii=1:height
    pix=fread(file_id,width,'float32');
   T13_real(ii,:)=pix;
end
fclose(file_id);

file_id=fopen('T3\T13_imag.bin','r');
status=fseek(file_id,0,'bof');
for ii=1:height
    pix=fread(file_id,width,'float32');
    T13_image(ii,:)=pix;
end
fclose(file_id);

file_id=fopen('T3\T33.bin','r');
status=fseek(file_id,0,'bof');
for ii=1:height
    pix=fread(file_id,width,'float32');
    T33(ii,:)=pix;
end
fclose(file_id);
% matlab文件操作结束

T13=complex(T13_real,T13_image);
Span = T11 + T22 + T33;
SpanMax = max(Span(:));
SpanMin = min(Span(:));
if SpanMin < eps
    SpanMin = eps;
end

Pv = zeros(height,width);
Ps = zeros(height,width);
Pd = zeros(height,width);

for ii = 1 : height
    for jj = 1 : width
        ALP = 0;
        BET = 0;
        TT22 = T22(ii,jj);
        TT13_im = T13_image(ii,jj);
        fv = 3 * TT22 / 2;
        TT11 = T11(ii,jj) - fv;
        TT33 = T33(ii,jj) - fv;
        TT13_re = T13_real(ii,jj) - fv / 3;
        if (TT11 <= eps) || (TT33 <= eps)      %Volume Scatter>Total
            fv = 3 * (TT11 + TT22 + TT33 + 2 *fv) /8 ;
            fd = 0;
            fs = 0;
        else
            %Data conditioning for non realizable Shh Svv* term
            if ((TT13_re * TT13_re + TT13_im * TT13_im) > TT11 * TT33)
                rtemp = TT13_re * TT13_re + TT13_im * TT13_im;
                TT13_re = TT13_re * sqrt(TT11 * TT33 / rtemp);
                TT13_im = TT13_im * sqrt(TT11 * TT33 / rtemp);
            end
            %Odd Bounce
            if (TT13_re >= 0)
                ALP = -1;
                fd = (TT11 * TT33 - TT13_re * TT13_re - TT13_im * TT13_im) / (TT11 + TT33 + 2 * TT13_re);
                fs = TT33 - fd;
                BET = sqrt((fd + TT13_re) * (fd + TT13_re) + TT13_im * TT13_im) / fs;
            %Even Bounce
            else
                BET = 1;
                fs = (TT11 * TT33 - TT13_re * TT13_re - TT13_im * TT13_im) / (TT11 + TT33 - 2 * TT13_re);
                fd = TT33 - fs;
                ALP = sqrt((fs - TT13_re) * (fs - TT13_re) + TT13_im * TT13_im) / fd;
            end
        end
        Pv(ii,jj) = 8 / 3 * fv;
        Pd(ii,jj) = fd *(1 + ALP^2);
        Ps(ii,jj) = fs *(1 + BET^2);
        
        Pv(ii,jj) = max(min(Pv(ii,jj),SpanMax),SpanMin);
        Pd(ii,jj) = max(min(Pd(ii,jj),SpanMax),SpanMin);
        Ps(ii,jj) = max(min(Ps(ii,jj),SpanMax),SpanMin);
    end
end

BLUE =imadjust(65535-uint16(-10*log(Ps)));
RED  =imadjust(65535-uint16(-10*log(Pd)));
GREEN=imadjust(65535-uint16(-10*log(Pv)));
RGB=cat(3,RED,GREEN,BLUE);

figure,imshow(BLUE,[149,150]),title('Ps');
saveas(1, 'Ps', 'tif');
figure,imshow(RED),title('Pd');
saveas(2, 'Pd', 'tif');
figure,imshow(GREEN),title('Pv');
saveas(3, 'Pv', 'tif');
figure,imshow(RGB),title('Freeman-Durden');
saveas(4, 'Freeman-decomposition', 'tif');

