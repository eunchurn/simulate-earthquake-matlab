close all;clear all; clc;

global hz z output_num data

%%%%% Data Processing %%%%
file1=char('wn_03.dat');

% Col.1 : Input Signal
% Col.2 : HMD Acc.
% Col.3 : 5th Top
% Col.4 : 5th Floor
% Col.5 : 2nd Floor
% Col.6 : 3rd Floor
% Col.7 : 4th Floor
% Col.8 : 2nd Disp.
g = 9.8; dt = 0.01;
nt = 2^13;
data=load(file1);
cal_coef1 = [1 g g   2 2 1 1 2.25/9/100]; % 비제어시

for ii=1:8
    data(:,ii)=data(:,ii)*cal_coef1(ii); % Calibration
end

for ii=1:6
[hz,frf(:,ii)]=exp_tf(data(:,1),data(:,ii+2),5,nt,dt); % transfer function
end

target=4;
tfrf=frf(:,6);

dt=0.01;
Y_elcen=load('elcacc1_analysis.dat');
Y_hachi=load('hachi_analysis.dat');
Y_kobe=load('kobe_analysis.dat');
Y_mexico=load('mexico_analysis.dat');
Y_north=load('north_analysis.dat');

Nstep=4096;%length(Y_elcen);
Nyf=1/(2*dt); 
if rem(Nstep,2)==0
    Nf=Nstep/2+1;
else
    Nf=(Nstep+1)/2;
end
df=1/(Nstep*dt);

YTarget=Y_elcen(:,1); % Select Target Response

Ys=fft(YTarget);
for ii=1:Nf
    if abs(tfrf(ii))<0.01*max(abs(tfrf))
        Us1(ii,1)=0;
    else
        Us1(ii,1)=Ys(ii,1)/tfrf(ii);
    end

end

if rem(Nstep,2)==0
   u=[Us1;conj(Us1(Nf-1:-1:2))];
else
  u=[Us1;conj(Us1(Nf:-1:2))];
end
t1=0:dt:(Nstep-1)*dt;
t2=0:dt:(length(Ys)-1)*dt;

u=ifft(u);u=real(u);
u=detrend(u);



load isg1.mat

% u=lsim(isg1,u,t1);
ec_elcen=u;
control_scale=5/max(ec_elcen)
figure('color',[1 1 1]),plot(t1,ec_elcen*control_scale)
save ec_elcen4.txt ec_elcen -ascii

edata=load('C:\Documents and Settings\Owner\바탕 화면\UNISON_0825\data\el_04.dat');
edata=dtrend(edata);
[bf,af,]=butter(7,10/50);
fedata=filtfilt(bf,af,edata);

figure('color',[1 1 1]),plot(t1,fedata(:,8)*5,t2,Y_elcen(:,1))
