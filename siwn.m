close all;clear all; clc;

global hz z output_num data

%%%%% Data Processing %%%%
file1=char('wn_01.dat');

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
cal_coef = [1 g g   2 2 1 1 2.25/9]; % 비제어시

for ii=1:8
    data(:,ii)=data(:,ii)*cal_coef(ii); % Calibration
end

for ii=1:6
    [hz,frf(:,ii)]=exp_tf(data(:,2),data(:,ii+2),5,nt,dt); % transfer function
end
lim_hz=[0.1 3];
lege=char(['5th top  ';'5th floor';'2nd floor';'3rd floor';'4th floor';'2nd disp.']);
figure('color',[1 1 1])
subplot(211),plot(hz,abs(frf(:,1)),':',hz,abs(frf(:,2)),'--',hz,abs(frf(:,3)),hz,abs(frf(:,4)),hz,abs(frf(:,5)),'LineWidth',2)
xlim(lim_hz),xlabel('Frequency'),ylabel('Magnitude'),legend(lege)
subplot(212),plot(hz,phase(frf(:,1)),':',hz,phase(frf(:,2)),'--',hz,phase(frf(:,3)),hz,phase(frf(:,4)),hz,phase(frf(:,5)),'LineWidth',2)

[hz_hmd,frf_hmd]=exp_tf(data(:,1),data(:,2),5,nt,dt); % transfer function
lim_hz2=[0 20];

uc_frf=[frf(:,3) frf(:,4) frf(:,5) frf(:,2) -frf(:,1)];

%%%% 모드추출 %%%

df2=(hz(2)-hz(1));
band=ceil([0.3 1.5 3]/df2); %FFT분석을 보고  대강의 주파수 대역
nmode=length(band)-1;
n=5;
for ii=1:nmode
    [a,b]=max(uc_frf(band(ii):band(ii+1),:),[],1);
    tf(:,ii)=a.';
    for jj=1:n
        fn1(jj,ii)=hz(band(ii)+b(jj)-1);
    end
end

fn2=mean(fn1);
TF=tf;
abs_tf=abs(TF);
angle_tf=angle(TF);
vangle=diff(angle_tf);

nmode=2;n=5;
exphi(1,:)=abs_tf(1,:);
for ii=1:nmode
    for jj=1:n-1
        if abs(vangle(jj,ii)) <pi/2
            exphi(jj+1,ii)=sign(exphi(jj,ii))*abs_tf(jj+1,ii);
        else
            exphi(jj+1,ii)=-sign(exphi(jj,ii))*abs_tf(jj+1,ii);
        end
    end
end

%%%%%%%%%%%%%% [1] Unison Tower의 약축방향의 모델(5DOF) %%%%%%%%%%%%%%%
%
% 1) Mass & Mass matrix
%
% Mass
m1_1 = 64800 / 5 ;                         % kg (층 베이스의 질량)
m2_1 = 32027.5 / 5 ;                      % kg (기둥질량+층 Frame)
m1 = m1_1 +  m2_1 ;                      % lumped mass1
m2 = m1_1 +  m2_1 ;                      % lumped mass2
m3 = m1_1 +  m2_1 ;                      % lumped mass3
m4 = m1_1 +  m2_1 ;                      % lumped mass4
m5 = m2_1 ;                             % lumped mass5
Total_mass = m1+m2+m3+m4+m5
% Mass matrix M
M = diag([m1 m2 m3 m4 m5]) ;
%
% 2) Modal matrix
%
% Mode shape getted from experiment
%
% from Ansys Analysis
%
MS(:,1) = [  0.16232      0.25691       0.33291         0.38769         0.41759  ]' ;   % 1st mode
MS(:,2) = [  0.38555      0.40494       0.17702       -0.15870        -0.40372  ]' ;   % 2nd mode
MS(:,3) = [  0.43733      0.030694    -0.42893      -0.20268          0.35121  ]' ;   % 3rd mode
MS(:,4) = [  0.36144    -0.39011      -0.074915      0.44445         -0.26250  ]' ;    % 4th mode
MS(:,5) = [  0.20417    -0.40316        0.46236      -0.36264          0.14184  ]' ;    % 5th mode
% Modal matrix W
tp = 0.0 ;
for i = 1:5
    for j = 1:5
        tp = tp + MS(j,i)^2 ;
    end
    W(:,i)  = MS(:,i)./sqrt(tp) ;
    tp = 0.0 ;
end

nMass=exphi'*M*exphi;
exphi=exphi*sqrt(diag(1./diag(nMass)))

nMass2=W'*M*W;
W=W*sqrt(diag(1./diag(nMass2)))

exMn=exphi'*M*exphi
% Gammak=exphi'*M*ones(n,1);

%% xi ID %%
md=1500;
B=[0;md;0;0;0];
Gammak=exphi'*B;
exwn=fn2*2*pi;
for ii=1:nmode
    for jj=1:n
        xik(jj,ii)=(2*sqrt((abs_tf(jj,ii)/exphi(jj,ii)/Gammak(ii))^2-1))^-1;
    end
end
figure('color',[1 1 1])
subplot(121),plot([0;exphi(:,1)],[0:5],[0;exphi(:,2)],[0:5]),legend('1st mode','2nd mode')
title('Measured model')

%   clear MS tp
%  Effective Mass
Mef = (W' * diag(M)).^2 ./ diag(W'*M*W);
fprintf(1,['\n Effective mass are \n\n',...
    '%7.1f kg, %7.1f kg, %7.1f kg, %7.1f kg, %7.1f kg %7.1f kg\n\n'],Mef);
% Total mass
Mtot = sum(diag(M))

%%% Estimate Stiffness matrix K & Damping matrix C
Z = diag(0.0198*[1  1  1  1  1]) ;                                         % from experiment(1% Damping)
fr  =  [0.5022,   1.5623,   2.74,   3.91,   4.83] ;                          % Hz(from Ansys)
fn  =  fr'./sqrt(1-diag(Z).^2) ;                                          % Hz
wn =  diag(2*pi*fn) ;

%K = W*(diag(diag([W'*M*W]))*(wn^2))*W'
%C = W*(diag(diag([W'*M*W]))*(2*Z*wn))*W'

K=inv(W')*(wn^2)*inv(W);
C=M*W*diag(diag(2*Z.*wn))*W'*M;

[Wn4,Phi4]=eign(K,M);
fn5=sqrt(Wn4)./(2*pi);

subplot(122),plot([0;Phi4(:,1)],[0:5],[0;Phi4(:,2)],[0:5]),legend('1st mode','2nd mode'),title('FE Updated model')

%%%% FE Model Update %%%%
%xik
%exK=M*exphi*diag(exwn.^2)*exphi'*M
%exC=M*exphi*diag(2*xik(1,:).*exwn)*exphi'*M

exLambda=diag(exwn)^2;

%% (Baruch and Barltzhack method) %%
upphi=exphi*(exphi'*M*exphi)^(-1/2);
upK1=K-(K*upphi*upphi'*M)-(M*upphi*upphi'*K)+(M*upphi*upphi'*K*upphi*upphi'*M)+(M*upphi*exLambda*upphi'*M);

[Wn2,Phi2]=eign(upK1,M);
fn3=sqrt(Wn2)./(2*pi)

%% (Berman and Nagy method) %%
Ma1=exphi'*M*exphi;
upM=M+M*exphi*inv(Ma1)*(eye(nmode)-Ma1)*inv(Ma1)*exphi'*M;
upK2=K-K*exphi*exphi'*M-M*exphi*exphi'*K+M*exphi*exphi'*K*exphi*exphi'*M+M*exphi*exLambda*exphi'*M;

[Phi3,Wn3]=eig(upK2,upM);
fn4=sqrt(Wn3)./(2*pi);

K2=upK1;   % Baruch et al. 선택

%%%% Modal Assurance Criterion %%%%
for ii=1:nmode
    for jj=1:nmode
        mac(ii,jj)=abs(exphi(:,jj)'*Phi2(:,ii))^2/((Phi2(:,ii)'*Phi2(:,ii))*(exphi(:,jj)'*exphi(:,jj)));
        mod_mac(ii,jj)=abs(exphi(:,jj)'*M*Phi2(:,ii))^2/((Phi2(:,ii)'*M*Phi2(:,ii))*(exphi(:,jj)'*M*exphi(:,jj)));
    end
end
disp('Modal Assurance Criterion (Orthogonality Check)')
mac
mod_mac
%%%% FE Model Updating 결과의 HMD 가속도입력을 받는 System 모델링 %%%

Aup=[zeros(5) eye(5);-inv(M)*K2 -inv(M)*C];
Bup=[zeros(5,1);inv(M)*B];
Cup=[-inv(M)*K2 -inv(M)*C];
Dup=inv(M)*B;

[numup,denup]=ss2tf(Aup,Bup,Cup,Dup);
fmax=8;
f = 0.001:0.001:fmax;
w = 2*pi*f;

h1 = freqs(numup(1,:), denup, w); mag1 = abs(h1); ph1 = angle(h1).*180/pi ;
h2 = freqs(numup(2,:), denup, w); mag2 = abs(h2); ph2 = angle(h2).*180/pi ;
h3 = freqs(numup(3,:), denup, w); mag3 = abs(h3); ph3 = angle(h3).*180/pi ;
h4 = freqs(numup(4,:), denup, w); mag4 = abs(h4); ph4 = angle(h4).*180/pi ;
h5 = freqs(numup(5,:), denup, w); mag5 = abs(h5); ph5 = angle(h5).*180/pi ;

figure('color',[1 1 1])
subplot(211),semilogy(hz,abs(frf(:,5)),':',f,mag5)
axis([0 fmax 0.0001 1])
title('Accelerance at the 5th floor')
subplot(212),semilogy(hz,abs(frf(:,4)),':',f,mag4)
axis([0 fmax 0.0001 1])
title('Accelerance at the 4th floor')

figure('color',[1 1 1])
subplot(311),semilogy(hz,abs(frf(:,3)),':',f,mag3)
axis([0 fmax 0.0001 1])
title('Accelerance at the 3rd floor')
subplot(312),semilogy(hz,abs(frf(:,2)),':',f,mag2)
axis([0 fmax 0.0001 1])
title('Accelerance at the 2nd floor')
subplot(313),semilogy(hz,abs(frf(:,1)),':',f,mag1)
axis([0 fmax 0.0001 1])
title('Accelerance at the 1st floor')

%%%%% Curve fitting HMD transfer function %%%%%
out=data(:,1);   % input / output data
in=data(:,2);

N = 10; re = 2^10;      % No. of reapetitive experiment and overlapped data
nt2=2^12;
NN1=5; shz=4; lhz=500;  % No. of order, first and last index for fitting SYS
code1=2;                % 0 : norestriction on pole location of SYS
                        % 1 : restrict SYS to be stable, minimum-phase
                        % 2 : constrain rational fit so that SYS is stable
NN2=11; code2=1;        % fitting of inverse SYS

T1 = [0:dt:dt*(N*nt2-1)]';
nf = nt/2; zf = 1:nf;                % number of frequency vector
dhz= 1/(nt*dt);                        % frequency interval
hz = [0:dhz:(nf-1)*dhz]'; hz(1)=0.001; % frequency vector
%%% HMD Transfer function curve fitting %%%
frdat= vpck(frf_hmd(zf),hz*2*pi);
frd  = xtract(frdat, hz(shz)*2*pi, hz(lhz)*2*pi);

sys  = fitsys(frd,NN1,1,code1,2,0);
[Ag1,Bg1,Cg1,Dg1] = unpck(sys);         % Transfer Function of the Excitor
sg1  = ss(Ag1,Bg1,Cg1,Dg1); sdg1 = c2d(sg1,dt,'zoh');
[Adg1,Bdg1,Cdg1,Ddg1]=ssdata(sdg1);       % Discrete model
[num_g1,den_g1]=ss2tf(Ag1,Bg1,Cg1,Dg1);
[Mag1, Phse1]=bode(num_g1,den_g1,hz*2*pi);
figure('color',[1 1 1])
subplot(2,1,1),semilogy(hz,abs(frf_hmd(zf)),'--',hz,Mag1,'-','LineWidth',2)
xlim(lim_hz2),xlabel('Frequency');ylabel('Magnitude');
title('HMD Transfer function system curve fitting result')
legend('Experimental TF','Curve Fitted TF')
subplot(2,1,2),plot(hz,phase(frf_hmd(zf))*(180/pi),'--',hz,Phse1+360,'-','LineWidth',2)
xlim(lim_hz2),xlabel('Frequency');ylabel('Phase');

%%% Inverse HMD Transfer function curve fitting %%%
ifrdat= vpck(1./frf_hmd(zf),hz*2*pi);
ifrd  = xtract(ifrdat, hz(shz)*2*pi, hz(lhz)*2*pi);

isys  = fitsys(ifrd,NN2,1,code2,2,0);
[iAg1,iBg1,iCg1,iDg1] = unpck(isys);         % Transfer Function of the Excitor
isg1  = ss(iAg1,iBg1,iCg1,iDg1); isdg1 = c2d(isg1,dt,'zoh');
[iAdg1,iBdg1,iCdg1,iDdg1]=ssdata(isdg1);       % Discrete model
[inum_g1,iden_g1]=ss2tf(iAg1,iBg1,iCg1,iDg1);
[Mag2, Phse2]=bode(inum_g1,iden_g1,hz*2*pi);
tc = lsim(iAg1,iBg1,iCg1,iDg1,in,T1);         % Simulated control signal of HMD
figure('color',[1 1 1])
subplot(2,1,1),semilogy(hz,abs(1./frf_hmd(zf)),'--',hz,Mag2,'-','LineWidth',2)
xlim(lim_hz2),xlabel('Frequency');ylabel('Magnitude');
title('Inverse HMD Transfer function system curve fitting result')
legend('Experimental TF','Curve Fitted TF')
subplot(2,1,2),plot(hz,phase(1./frf_hmd(zf))*(180/pi),'--',hz,Phse2+180,'-','LineWidth',2)
xlim(lim_hz2),xlabel('Frequency');ylabel('Phase');

%%%% Earthquake data Processing %%%%

neq=5; % Number of eq
scale=0.05;
eq1=load('eq\elcacc1.txt')./100*scale;  % dt=0.01 | cm/sec^2 -> m/sec^2
eq2=load('eq\hachi.txt')./100*scale;
eq3=load('eq\kobe.txt')./100*scale;
eq4=load('eq\mexico.txt')./100*scale;
eq5=load('eq\north.txt')./100*scale;
for ii=1:neq
    tii=num2str(ii);
    eval(['Nst(',tii,')=length(eq',tii,');'])
    eval(['t',tii,'=[0:dt:(Nst(',tii,')-1)*dt];'])
end
Nyf=1/(2*dt);
for ii=1:neq
    tii=num2str(ii);
    if rem(Nst(ii),2)==0
        Nf(ii)=Nst(ii)/2+1;
    else
        Nf(ii)=(Nst(ii)+1)/2;
    end
    df(ii)=1/(Nst(ii)*dt);
    eval(['Freq',tii,'=transpose([0:df(',tii,'):df(',tii,')*(Nf(',tii,')-1)]);'])
    eval(['Omega',tii,'=2*pi*Freq',tii,';'])
end

%%%% Earthquake Excitation %%%%
H=-M*ones(n,1);       % Base Excitation
Bxg=[zeros(n,1); inv(M)*H];
[Adxg,Bdxg]=c2d(Aup,Bxg,dt);
for jj=1:neq
    tii=num2str(jj);
    eval(['Xxg',tii,'=zeros(2*n,Nst(',tii,'));'])
    eval(['for ii=1:Nst(',tii,')-1 ','Xxg',tii,'(:,ii+1)=Adxg*Xxg',tii,'(:,ii)+Bdxg*eq',tii,'(ii);','end'])
    eval(['Yxg',tii,'=Cup*Xxg',tii,';'])
end

%%%% 지반 가진 Target 층의 가속도를 실현하는 HMD입력신호 생성
Tar_Resp=5;
YTarget1=Yxg1(Tar_Resp,:)'; % El Centro Response
YTarget2=Yxg2(Tar_Resp,:)'; % Hachinohe Response
YTarget3=Yxg3(Tar_Resp,:)'; % Kobe Response
YTarget4=Yxg4(Tar_Resp,:)'; % Mexico City Response
YTarget5=Yxg5(Tar_Resp,:)'; % Northridge Response
Cxg=Aup(n+Tar_Resp,:);Dxg=zeros(1,1);

% El Centro
Cupt=Cup(Tar_Resp,:);Dupt=Dup(Tar_Resp,:);
[num,dem]=ss2tf(Aup,Bup,Cupt,Dupt);
[Tyu]=freqs(num,dem,Freq1*2*pi);
Ys1=fft(YTarget1);
for ii=1:Nf(1)
    if abs(Tyu(ii,1))<0.05*max(abs(Tyu))
        Us1(ii,1)=0;
    else
        Us1(ii,1)=Ys1(ii,1)/Tyu(ii,1);
    end
end

if rem(Nst(1),2)==0
    u1=[Us1;conj(Us1(Nf(1)-1:-1:2))];
else
    u1=[Us1;conj(Us1(Nf(1):-1:2))];
end
aTyu=abs(Tyu);

u1=ifft(u1);u1=real(u1);
u1=detrend(u1);

con_eq1=lsim(iAg1,iBg1,iCg1,iDg1,u1,t1); % Generate HMD control signal
ret_resp=lsim(Ag1,Bg1,Cg1,Dg1,con_eq1,t1); % HMD response

Yhmd1=lsim(Aup,Bup,Cup,Dup,u1,t1);  % HMD Excited
figure('color',[1 1 1])
for ii=1:5
    subplot(5,1,ii),plot(t1,Yxg1(ii,:),t1,Yhmd1(:,ii))
    title([num2str(ii),' Story'])
end
legend('Earthquake Excited','HMD Excited')

figure('color',[1 1 1])
subplot(211),plot(t1,con_eq1),title('control signal of HMD')
subplot(212),plot(t1,u1,t1,-ret_resp),title('Estimate HMD Response'),legend('Exact','HMD Response')

save control_signals\cont_elcen.txt con_eq1 -ascii
