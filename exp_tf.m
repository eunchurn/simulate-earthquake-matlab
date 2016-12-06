% exp_tf - Experimental transfer function
%
% USAGE ::
% Return Frequency Response Function of Data file
%
% [frf,hz] = exp_tf(data,inc,ouc,N,nt,dt,cal_coef)  
%
% in_data = input column data
% out_data = output column data
% inc = input channel (column vector)
% ouc = output channel (column vector)
% N = Reapetitive number of experiment 
% nt = Number of Time Vector
% dt = Time intercal


function [hz,frf] = exp_tf(in_data,out_data,N,nt,dt) 

global M T in Exp Ana sector dhz z output_num data nt 

 re = nt/2;   % averaging
 
%%%====================%%%
%%%= Experimental FRF =%%%
%%%====================%%%  

 [nl,rl]  = size(in_data);     %total data length
 zl  = nl-(nt*N)+1:nl;   %used data length 
 in_ran=dtrend(in_data(zl)); 
 out_ran=dtrend(out_data(zl)); 
 
 t=0:dt:dt*(nt-1);    %time vector
 t1=0:dt:dt*(N*nt-1);
 nf=nt/2; z=1:nf;     %number of frequency vector
 dhz=1/(nt*dt);       %frequency interval
 hz=0:dhz:(nt-1)*dhz; %frequency vector

 in=in_ran(1:nt*N);  % input data sequency
 out=out_ran(1:nt*N); % output data sequency

% Obtain the transfer function 
 I=0;  
 fft_in=zeros(nt,1); 
 fft_out=zeros(nt,1);
 frf_in=zeros(nt,1); 
 frf_out=zeros(nt,1); 
 frf=zeros(nt,1);
 while re*(I-1)+nt < nt*N
    I=I+1;
    
    fft_in  = fft( in(re*(I-1)+1 : re*(I-1)+nt , 1)); 
    % FFT the input data  (input voltage)
    fft_out = fft(out(re*(I-1)+1 : re*(I-1)+nt , 1)); 
    % FFT the output data (output acceleration)
       
    frf_in  = frf_in + fft_in .* conj(fft_in);
    frf_out = frf_out + fft_out.* conj(fft_in);
 end
 frf = frf_out ./ frf_in;
