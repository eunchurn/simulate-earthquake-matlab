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
cal_coef1 = [1 g g   2 2 1 1 2.25/9]; % 비제어시
cal_coef2 = [1 g g   2 2 1 1 2.25/9]; % 비제어시
for ii=1:8
    data(:,ii)=data(:,ii)*cal_coef1(ii); % Calibration
end

for ii=1:6
    [hz,frf(:,ii)]=exp_tf(data(:,2),data(:,ii+2),5,nt,dt); % transfer function
end
lim_hz=[0.3 8];
lege=char(['5th top  ';'5th floor';'2nd floor';'3rd floor';'4th floor';'2nd disp.']);
figure('color',[1 1 1])
subplot(211),plot(hz,abs(frf(:,1)),':',hz,abs(frf(:,2)),'--',hz,abs(frf(:,3)),hz,abs(frf(:,4)),hz,abs(frf(:,5)),'LineWidth',2)
xlim(lim_hz),xlabel('Frequency'),ylabel('Magnitude'),legend(lege)
subplot(212),plot(hz,phase(frf(:,1))*(180/pi)+360,':',hz,phase(frf(:,2))*(180/pi)+210,'--',hz,phase(frf(:,3))*(180/pi)+900,hz,phase(frf(:,4))*(180/pi)-180,hz,phase(frf(:,5))*(180/pi)+900,'LineWidth',2)
xlim(lim_hz),xlabel('Frequency'),ylabel('Phase (degree)')