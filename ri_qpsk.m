clc
clear all
close all

%BPSK
A=5;
t=0:.001:1;
f1=10;
f2=2;
b=A.*sin(2*pi*f1*t);%Carrier Sine
subplot(3,1,1);
%figure
plot(t,b);
xlabel('time');
ylabel('Amplitude');
title('Carrier');
grid on;

u=square(2*pi*f2*t);%Message signal
subplot(3,1,2);
%figure
plot(t,u);
xlabel('time');
ylabel('Amplitude');
title('Message Signal');
grid on;

v=b.*u;
subplot(3,1,3);
plot(t,v)

%QPSK
data=[0  1 0 1 1 1 0 0 1 1]; % information

data_NZR=2*data-1; % Data Represented at NZR form for QPSK modulation
s_p_data=reshape(data_NZR,2,length(data)/2);  % S/P convertion of data

br=10.^6; %transmission bit rate  
f=br;
T=1/br; % bit duration
t=T/99:T/99:T; % Time vector for one bit information

%QPSK modulation
y=[];
y_in=[];
y_qd=[];
for(o=1:length(data)/2)
    y1=s_p_data(1,o)*cos(2*pi*f*t); % inphase component
    y2=s_p_data(2,o)*sin(2*pi*f*t) ;% Quadrature component
    y_in=[y_in y1]; % inphase signal vector
    y_qd=[y_qd y2]; %quadrature signal vector
    y=[y y1+y2]; % modulated signal vector
end
Tx_sig=y; % transmitting signal after modulation
tt=T/99:T/99:(T*length(data))/2;

figure(2)

subplot(3,1,1);
plot(tt,y_in,'linewidth',3), grid on;
title(' wave form for inphase component in QPSK modulation ');
xlabel('time(sec)');
ylabel(' amplitude(volt0');

subplot(3,1,2);
plot(tt,y_qd,'linewidth',3), grid on;
title(' wave form for Quadrature component in QPSK modulation ');
xlabel('time(sec)');
ylabel(' amplitude(volt0');

subplot(3,1,3);
plot(tt,Tx_sig,'r','linewidth',3), grid on;
title('QPSK modulated signal (sum of inphase and Quadrature phase signal)');
xlabel('time(sec)');
ylabel(' amplitude(volt0');



m=2;
L = 2*m;%Nyquist equation

T=0.01:0.001:.05;
Pf =0.01:0.01:1;
snr_dB = -4;
snr_dB1 = 0;
snr_dB2 = 2;
%-----Linear Value of SNR-----%
snr1 = 10.^(snr_dB./20);
snr2 = 10.^(snr_dB1./20);
snr3 = 10.^(snr_dB2./20);
for p = 1:length(Pf)
    g = 0;
    j = 0;
    k = 0;
    l = 0;
    c = 0;
    d = 0;
       
    for kk=1:400 % No. of Monte Carlo Simulations
        Noise1=1/sqrt(2)*(randn(1,L)+1i*randn(1,L));  % AWGN noise with mean=0 var=1
        Noise2=1/sqrt(2)*(randn(1,L)+2i*randn(1,L));
        Noise3=1/sqrt(2)*(randn(1,L)+3i*randn(1,L));
        
       K=2; %Rician fading
         mu = sqrt( K/(2*(K+1))); 
         su = sqrt( 1/(2*(K+1)));
         ri = ( su*randn(1,L) + mu ) + 1i*( su*randn(1,L) + mu );
         R = abs(ri);
         
       
        
        n1=Noise1*snr1;
        n2=Noise2*snr2;
        n3=Noise3*snr3;
        
        pn1 = (std(Noise1)).^2;%Noise power 
        pn2 = (std(Noise2)).^2;
        pn3 = (std(Noise3)).^2;
        
        Signal = v(1:L); %BPSK
        Recv_Sig1 = R.*Signal + n1;
        Recv_Sig2 = R.*Signal + n2;
        Recv_Sig3 = R.*Signal + n3;
        
        Sign = Tx_sig(1:L);%QPSK
        Rev_Sig1 = R.*Sign + n1;
        Rev_Sig2 = R.*Sign + n2;
        Rev_Sig3 = R.*Sign + n3;
        
        % Energy of received signal
        Energy1 = abs(Recv_Sig1).^2; %BPSK
        Energy2 = abs(Recv_Sig2).^2;
        Energy3 = abs(Recv_Sig3).^2;
        
        Energy4 = abs(Rev_Sig1).^2; %QPSK
        Energy5 = abs(Rev_Sig2).^2;
        Energy6 = abs(Rev_Sig3).^2;
        
        %Threshold Calculation
        Threshold(p) = gaminv(1-Pf(p),m)*2;
         
        
       
        signal_power1(p) = sum(abs(Recv_Sig1.^2))/pn1; %BPSK
        signal_power2(p) = sum(abs(Recv_Sig2.^2))/pn2;
        signal_power3(p) = sum(abs(Recv_Sig3.^2))/pn3;
        signal_power4(p) = sum(abs(Rev_Sig1.^2))/pn1; %QPSK
        signal_power5(p) = sum(abs(Rev_Sig2.^2))/pn2;
        signal_power6(p) = sum(abs(Rev_Sig3.^2))/pn3;
        
        
        if(signal_power1(p) > Threshold(p))  %BPSK
            j=j+1;
        end
        if(signal_power4(p) > Threshold(p))  %BPSK
            g=g+1;
        end
        
        if(signal_power2(p) > Threshold(p))
            l=l+1;
        end
        if(signal_power5(p) > Threshold(p))
            k=k+1;
        end
        
        if(signal_power3(p) > Threshold(p))  
            d=d+1;
        end
        if(signal_power6(p) > Threshold(p))  
            c=c+1;
        end
end
Pd1(p) = j/kk;
Pd2(p)= g/kk;
Pd3(p) = l/kk;
Pd4(p) = k/kk;
Pd5(p) = d/kk;
Pd6(p) = c/kk;

pm1(p)=1-Pd1(p);
pm2(p)=1-Pd2(p);
pm3(p)=1-Pd3(p);
pm4(p)=1-Pd4(p);
pm5(p)=1-Pd5(p);
pm6(p)=1-Pd6(p);


end
figure
loglog(Pf,pm1,'b',Pf,pm3,'r',Pf,pm5,'c','LineWidth', 3);
hold on
grid on
xlabel('Probability of false alarm,P_{f}','FontSize',30);
ylabel('Probability of misdetection,P_{m}','FontSize',30);
title('CROC Curve Under Rician Channel','FontSize',30);
leg=legend('Rician for snr= -4dB','Rician for snr= 0dB', 'Rician for snr= 2dB');
set(leg,'FontSize',30)