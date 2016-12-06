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
cal_coef1 = [1 g g   2 2 1 1 2.25/9]; % 비제어시

for ii=1:8
    data(:,ii)=data(:,ii)*cal_coef1(ii); % Calibration
end

for ii=1:6
[hz,frf(:,ii)]=exp_tf(data(:,1),data(:,ii+2),5,nt,dt); % transfer function
end
%%%%% Curve fitting HMD transfer function %%%%%
out=data(:,1);   % input / output data
in=data(:,2);

N = 10; re = 2^12;      % No. of reapetitive experiment and overlapped data
nt2=2^12;
NN1=15; shz=4; lhz=500;  % No. of order, first and last index for fitting SYS
code1=1;                % 0 : norestriction on pole location of SYS
                        % 1 : restrict SYS to be stable, minimum-phase
                        % 2 : constrain rational fit so that SYS is stable
NN2=11; code2=1;        % fitting of inverse SYS
lim_hz2=[0 20];

T1 = [0:dt:dt*(N*nt2-1)]';
nf = nt/2; zf = 1:nf;                % number of frequency vector
dhz= 1/(nt*dt);                        % frequency interval
hz = [0:dhz:(nf-1)*dhz]'; hz(1)=0.001; % frequency vector
%%% HMD Transfer function curve fitting %%%
frdat= vpck(frf(zf,1),hz*2*pi);
frd  = xtract(frdat, hz(shz)*2*pi, hz(lhz)*2*pi);

sys  = fitsys(frd,NN1,1,code1,2,0);
[Ag1,Bg1,Cg1,Dg1] = unpck(sys);         % Transfer Function of the Excitor
sg1  = ss(Ag1,Bg1,Cg1,Dg1); sdg1 = c2d(sg1,dt,'zoh');
[Adg1,Bdg1,Cdg1,Ddg1]=ssdata(sdg1);       % Discrete model
[num_g1,den_g1]=ss2tf(Ag1,Bg1,Cg1,Dg1);
[Mag1, Phse1]=bode(num_g1,den_g1,hz*2*pi);
figure('color',[1 1 1])
subplot(2,1,1),semilogy(hz,abs(frf(zf,1)),'--',hz,Mag1,'-','LineWidth',2)
xlim(lim_hz2),xlabel('Frequency');ylabel('Magnitude');
title('HMD Transfer function system curve fitting result')
legend('Experimental TF','Curve Fitted TF')
subplot(2,1,2),plot(hz,phase(frf(zf,1))*(180/pi),'--',hz,Phse1+360,'-','LineWidth',2)
xlim(lim_hz2),xlabel('Frequency');ylabel('Phase');

%%% Inverse HMD Transfer function curve fitting %%%
ifrdat= vpck(1./frf(zf,1),hz*2*pi);
ifrd  = xtract(ifrdat, hz(shz)*2*pi, hz(lhz)*2*pi);

isys  = fitsys(ifrd,NN2,1,code2,2,0);
[iAg1,iBg1,iCg1,iDg1] = unpck(isys);         % Transfer Function of the Excitor
isg1  = ss(iAg1,iBg1,iCg1,iDg1); isdg1 = c2d(isg1,dt,'zoh');
[iAdg1,iBdg1,iCdg1,iDdg1]=ssdata(isdg1);       % Discrete model
[inum_g1,iden_g1]=ss2tf(iAg1,iBg1,iCg1,iDg1);
[Mag2, Phse2]=bode(inum_g1,iden_g1,hz*2*pi);
[Mag3, Phse3]=bode(den_g1,num_g1,hz*2*pi);
tc1 = lsim(iAg1,iBg1,iCg1,iDg1,in,T1);         % Simulated control signal of HMD
[iAg2,iBg2,iCg2,iDg2]=tf2ss(den_g1,num_g1);
tc2 = lsim(iAg2,iBg2,iCg2,iDg2,in,T1);

figure('color',[1 1 1])
subplot(211),plot(T1,out,T1,tc1)
subplot(212),plot(T1,out,T1,tc2),xlim([0 10])

figure('color',[1 1 1])
subplot(2,1,1),semilogy(hz,abs(1./frf(zf,1)),'--',hz,Mag2,'-','LineWidth',2)
xlim(lim_hz2),xlabel('Frequency');ylabel('Magnitude');
title('Inverse HMD Transfer function system curve fitting result')
legend('Experimental TF','Curve Fitted TF')
subplot(2,1,2),plot(hz,phase(1./frf(zf,1))*(180/pi),'--',hz,Phse2+180,'-','LineWidth',2)
xlim(lim_hz2),xlabel('Frequency');ylabel('Phase');

figure('color',[1 1 1])
subplot(2,1,1),semilogy(hz,abs(1./frf(zf,1)),'--',hz,Mag3,'-','LineWidth',2)
xlim(lim_hz2),xlabel('Frequency');ylabel('Magnitude');
title('Inverse HMD Transfer function system curve fitting result')
legend('Experimental TF','Curve Fitted TF')
subplot(2,1,2),plot(hz,phase(1./frf(zf,1))*(180/pi),'--',hz,Phse3+180,'-','LineWidth',2)
xlim(lim_hz2),xlabel('Frequency');ylabel('Phase');

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

YTarget=Y_elcen(:,15); % Select Target Response

Ys=fft(YTarget);
for ii=1:Nf
    if abs(frf(ii,1))<0.01*max(abs(frf(:,1)))
        Us1(ii,1)=0;
    else
        Us1(ii,1)=Ys(ii,1)/frf(ii,1);
    end

end

if rem(Nstep,2)==0
   u=[Us1;conj(Us1(Nf-1:-1:2))];
else
  u=[Us1;conj(Us1(Nf:-1:2))];
end
t1=0:dt:(Nstep-1)*dt;
t2=0:dt:(length(Ys)-1)*dt;
edata=load('C:\Documents and Settings\Owner\바탕 화면\UNISON_0825\data\el_02.dat');
edata=dtrend(edata);
[bf,af,]=butter(7,10/50);
fedata=filtfilt(bf,af,edata);

plot(t1,fedata(:,3)*5,t2,Y_elcen(:,15))

u=ifft(u);u=real(u);
u=detrend(u);

load isg1.mat

% u=lsim(isg1,u,t1);
ec_elcen=u;
save ec_elcen4.txt ec_elcen -ascii
