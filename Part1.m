% Ryan Lindsay 
% 101038101

%% Part 1

% Most of the code is taken directly from the MNPA 
% Givens
R1 = 1;
R2 = 2;
R4 = 0.1;
RO = 1000;


R3 = 20;

G1 = 1/R1;
G2 = 1/R2;
G3 = 1/R3;
G4 = 1/R4;
GO = 1/RO;

L = 0.2;
Cap = 0.25;
alpha = 100;
 
%Matrix Definitions

C = [0 0 0 0 0 0 0;
   -Cap Cap 0 0 0 0 0;
    0 0 -L 0 0 0 0;
    0 0 0 0 0 0 0;
    0 0 0 0 0 0 0;
    0 0 0 0 0 0 0;
    0 0 0 0 0 0 0;];

G = [1 0 0 0 0 0 0;
   -G2 G1+G2 -1 0 0 0 0;
    0 1 0 -1 0 0 0;
    0 0 -1 G3 0 0 0;
    0 0 0 0 -alpha 1 0;
    0 0 0 G3 -1 0 0;
    0 0 0 0 0 -G4 G4+GO];


VDC=zeros(7,1);
VAC=zeros(7,1);
F=zeros(7,1);



for v = -10:0.1:10
    F(1,1)=v;
    VDC=G\F;
    
    figure(1)
    plot (v, VDC(7,1),'b.') 
    plot(v, VDC(4,1), 'g.') 
    hold on
end
    title('DC Sweep')
    xlabel('Vin (V)')
    ylabel('Voltage (V)')
    legend('V3','VO')
    
    hold off

w = logspace(0,5,500);                  
F(1,1) = 1;

for i = 1:length(w)
    VAC = (G+C*1j*w(i))\F;              
    figure(3)
    semilogx(w(i), abs(VAC(7,1)), 'b.')
    hold on
end

    xlabel('log (w)')
    ylabel('VO (V)')
    title('AC Sweep')
    hold off

for i = 1:length(w)
    VAC = (G+C*1j*w(i))\F; 
    gain = 20*log(abs(VAC(7,1))/F(1));   
    
    figure(4)
    plot(i, gain, 'g.')
    hold on 
end
    title('Gain Vo/Vin ')
    xlabel('Step')
    ylabel('Gain (dB)')
    hold off



perb =  Cap + 0.05.*randn(1,1000);
w = pi;
Gain = zeros(1000,1);

for n = 1:length(Gain)
    C(1,1) = perb(n);
    C(1,2) = -perb(n);
    C(2,1) = -perb(n);
    C(2,2) = perb(n);
    VAC = (G+C*1j*w)\F;                 
    Gain(n,1) = abs(VAC(7,1))/F(1);    
end


% figure(5)
% histogram(Gain,100);
% title('Histogram of Gain ')
% xlabel('Gain')
% ylabel('Number of elements')


%% Transient Circuit Simulation

Vin = 1;
steps = 1000;
dt = 1/steps;




F = [Vin; 0; 0; 0; 0; 0; 0;];

V0 = zeros(7,steps);
V1 = zeros(7,1);

Xval = 1:steps;

% Sim A

for i = 1:steps
    
    if i <30
        
        V0(:,i) = 0;
    
    elseif i == 30
        
        V0(:,i) = (C./dt+G)\(F+C*V1/dt);
        
    else
        
        V0(:,i) = (C./dt+G)\(F+C*Vprev/dt); 
        
    end
    
    Vprev = V0(:,i);
         
end

figure(6) 
plot(Xval, V0(7,:),'b')
hold on 
plot(Xval, V0(1,:),'r')
title(' Plot of Vin and Vout')
xlabel('Time (ms)')
ylabel('Voltage')
legend('Vout','Vin')



% Sim B 

V2 = zeros(7,steps);
freq = zeros(7,1);


Vprev = zeros(7,1);

for i =1:steps

freq(1) = sin(2*pi*(1/0.03)*i/steps);
V2(:,i) = (C./dt+G)\(freq+C*Vprev/dt);

Vprev = V2(:,i);

end

figure(7)
plot(Xval, V2(7,:), 'b')
hold on
plot(Xval, V2(1,:), 'r')
title('Vin vs Vout (Sine Pulse)')
xlabel('Time (ms)')
ylabel('Voltage')
legend('Vout','Vin')


% Sim C

V3 = zeros(7,steps);

Gpulse = zeros(7,1);


Vprev = zeros(7,1);

for i = 1:steps
    
    
    Gpulse(1,1) = exp(-1/2*((i/steps-0.06)/(0.03))^2);
    
    V3(:,i) = (C./dt+G)\(Gpulse+C*Vprev/dt);
    
    Vprev = V3(:,i);
    
end



figure(8)
plot(Xval, V3(7,:), 'b')
hold on
plot(Xval, V3(1,:), 'r')
title(' Vin vs Vout (Gaussian Pulse)')
xlabel(' Time (ms)')
ylabel('Voltage')
legend('Vout','Vin')
    
    
%% Fourier Transform

% Sim A
freq = linspace(-500,500,steps);

FV1 = fft(V0.');
FSV1 = fftshift(FV1);

figure(9)
plot(freq, 20*log10(abs(FSV1(:,1))),'r')
hold on 
plot(freq, 20*log10(abs(FSV1(:,7))),'b')
axis([-500 500 -100 inf])
ylabel('Voltage')
xlabel('Frequency')
title('Fourier Transform Plot')

% Sim B

FV2 = fft(V2.');
FSV2 = fftshift(FV2);

figure(10)
plot(freq, 20*log10(abs(FSV2(:,1))),'r')
hold on 
plot(freq, 20*log10(abs(FSV2(:,7))),'b')
axis([-500 500 -100 inf])
ylabel('Voltage')
xlabel('Frequency')
title('Fourier Transform Plot (Sine Pulse)')
% Sim C 

FV3 = fft(V3.');
FSV3 = fftshift(FV3);

figure(11)
plot(freq, 20*log10(abs(FSV3(:,1))),'r')
hold on 
plot(freq, 20*log10(abs(FSV3(:,7))),'b')
axis([-500 500 -100 inf])
ylabel('Voltage')
xlabel('Frequency')
title('Fourier Transform Plot (Gaussian Pulse')









