clc
clear all
close all
A=5;
t=0:.001:1;
f1=10;
f2=2;
y=A.*sin(2*pi*f1*t);%Carrier Sine
subplot(3,1,1);
%figure
plot(t,y);
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

v=y.*u;
subplot(3,1,3);
plot(t,v)

z=2.5;
L = 2*z;%Nyquist equation

t=0.01:0.001:.05;
Pf =0.01:0.01:1;
snr_dB =-6;
snr_dB1 =0;
snr_dB2 = 2;
%-----Linear Value of SNR-----%
snr1 = 10.^(snr_dB./10);
snr2 = 10.^(snr_dB1./10);
snr3 = 10.^(snr_dB2./10);
tic;
t1=tic;
for p = 1:length(Pf)
    g = 0;
    j = 0;
    k = 0;
    l = 0;
    c = 0;
    d = 0;
    m = 0;
    n = 0;
    o = 0;
    a=0;
    b=0;
    f=0;
       
    for kk=1:80000 % No. of Monte Carlo Simulations
        %x=sin(2*pi*f*t);
        Noise1=1/sqrt(2)*(randn(1,L)+1i*randn(1,L));  % AWGN noise with mean=0 var=1
        Noise2=1/sqrt(2)*(randn(1,L)+2i*randn(1,L));
        Noise3=1/sqrt(2)*(randn(1,L)+3i*randn(1,L));
        
        %rayleigh
        h1=1/sqrt(2)*(randn(1,L)+4i*randn(1,L));
        h2=1/sqrt(2)*(randn(1,L)+5i*randn(1,L));
        h3=1/sqrt(2)*(randn(1,L)+6i*randn(1,L));
        
        K=2; %Rician fading
         mu = sqrt( K/(2*(K+1))); 
         su = sqrt( 1/(2*(K+1)));
         ri = ( su*randn(1,L) + mu ) + 1i*( su*randn(1,L) + mu );
         R = abs(ri);
         
         %K=1; %Nakagami fading
         %s=sqrt(1/(2*(K+1)));
         %na=(s*randn(1,L)) + 1i*(s*randn(1,L));
         %N=abs(na);
      
        n1=Noise1*snr1;
        n2=Noise2*snr2;
        n3=Noise3*snr3;
        
        pn1 = (std(Noise1)).^2;%Noise power 
        pn2 = (std(Noise2)).^2;
        pn3 = (std(Noise3)).^2;
       
        Signal = v(1:L);
        Rev_Sig1 = h1.*Signal + n1; %rayleigh
        Rev_Sig2 = h2.*Signal + n2;
        Rev_Sig3 = h3.*Signal + n3;
        
        Recv_Sig1 = R.*Signal + n1; %rician
        Recv_Sig2 = R.*Signal + n2;
        Recv_Sig3 = R.*Signal + n3;
        
       
        
        
        Energy1 = abs(Rev_Sig1).^2; % Energy of received signal
        Energy2 = abs(Rev_Sig2).^2;
        Energy3 = abs(Rev_Sig3).^2;
        Energy4 = abs(Recv_Sig1).^2;
        Energy5 = abs(Recv_Sig2).^2;
        Energy6 = abs(Recv_Sig3).^2;
        
        
        %Threshold Calculation
        Threshold(p) = gaminv(1-Pf(p),z)*2;
        
       Signal_power1(p) = sum(abs(Rev_Sig1.^2))/pn1;
        signal_power1(p) = sum(abs(Recv_Sig1.^2))/pn1;
        
        Signal_power2(p) = sum(abs(Rev_Sig2.^2))/pn2;
        signal_power2(p) = sum(abs(Recv_Sig2.^2))/pn2;
        
        Signal_power3(p) = sum(abs(Rev_Sig3.^2))/pn3;
        signal_power3(p) = sum(abs(Recv_Sig3.^2))/pn3;
        
        if(Signal_power1(p) > Threshold(p))  
            g=g+1;
        end
        if(signal_power1(p) > Threshold(p))  
            j=j+1;
        end
        
        
        if(Signal_power2(p) > Threshold(p))  
           k=k+1;
        end
        if(signal_power2(p) > Threshold(p))  
            l=l+1;
        end
        
        if(Signal_power3(p) > Threshold(p))  
            c=c+1;
        end
        if(signal_power3(p) > Threshold(p))  
            d=d+1;
        end
        
end
Pd1(p) = g/kk;
Pd2(p)= j/kk;

Pd5(p) = k/kk;
Pd6(p) = l/kk;

Pd9(p) = c/kk;
Pd10(p) = d/kk;

pm1(p)= 1-Pd1(p);
pm2(p)= 1-Pd2(p);
pm5(p)= 1-Pd5(p);
pm6(p)= 1-Pd6(p);
pm9(p)= 1-Pd9(p);
pm10(p)= 1-Pd10(p);


end
t2=toc(t1);
t1;
t2;
figure
semilogy(Pf, Pd1,'-r',Pf,Pd2,'-*r',Pf,Pd5,'-g',Pf,Pd6,'-*g',Pf,Pd9,'-m',Pf,Pd10,'-*m','LineWidth', 3);
hold on
grid on
xlabel('Probability of false alarm,Pf','FontSize',12);
ylabel('Probability of detection,Pd','FontSize',12);
title('ROC CURVE UNDER Rayleigh and Rician  Channel','FontSize',12);
leg=legend('rayleigh for snr= -6db','rician for snr = -6db', 'rayleigh for snr= 0db ', 'rician for snr= 0 db', 'rayleigh for snr= 2db','rician for snr= 2db');
set(leg,'FontSize',12)