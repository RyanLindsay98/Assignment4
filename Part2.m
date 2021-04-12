%% Part 2: Circuit with Noise 

% Ryan Lindsay 
% 101038101
 
%Givens
R1 = 1;
R2 = 2;
R4 = 0.1;
RO = 1000;

% From Assignment 3

R3 = 20;


G1 = 1/R1;
G2 = 1/R2;
G3 = 1/R3;
G4 = 1/R4;
GO = 1/RO;

L = 0.2;
Cap = 0.25;
alpha = 100;


%Capacitor Values
Cn1 = 0.00001;

Cn2 = Cn1*2;

Cn3 = Cn1/1000;
%Cn3 = 0;


%Matrix Definitions

C1 = [0 0 0 0 0 0 0;
    -Cap Cap 0 0 0 0 0;
     0 0 -L 0 0 0 0;
     0 0 0 -Cn1 0 0 0;
     0 0 0 0 0 0 0;
     0 0 0 -Cn1 0 0 0;
     0 0 0 0 0 0 0;];

C2 = [0 0 0 0 0 0 0;
    -Cap Cap 0 0 0 0 0;
     0 0 -L 0 0 0 0;
     0 0 0 -Cn2 0 0 0;
     0 0 0 0 0 0 0;
     0 0 0 -Cn2 0 0 0;
     0 0 0 0 0 0 0;];

C3 = [0 0 0 0 0 0 0;
    -Cap Cap 0 0 0 0 0;
     0 0 -L 0 0 0 0;
     0 0 0 -Cn3 0 0 0;
     0 0 0 0 0 0 0;
     0 0 0 -Cn3 0 0 0;
     0 0 0 0 0 0 0;];

G = [1 0 0 0 0 0 0;
   -G2 G1+G2 -1 0 0 0 0;
    0 1 0 -1 0 0 0;
    0 0 -1 G3 0 0 0;
    0 0 0 0 -alpha 1 0;
    0 0 0 G3 -1 0 0;
    0 0 0 0 0 -G4 G4+GO];

Vin = 1;

F = [Vin; 0; 0; 0; 0; 0; 0;];


steps = 1000;

Xval = 1:steps;


V1 = zeros(7,1);
Gpulse = zeros(7,1);
Vprev = zeros(7,1);

dt = 1/steps;





% Sim A

for i = 1:steps
    
    % Current Source
    In = 0.001*randn();
    
    
    Gpulse(4,1) = In;
    
    Gpulse(1,1) = exp(-1/2*((i/steps-0.06)/(0.03))^2);
         
    V1(:,i) = (C1./dt+G)\(Gpulse+C1*Vprev/dt);
    
    Vprev = V1(:,i);
    
end

figure(1)
plot(Xval, V1(7,:),'b')
hold on
plot(Xval, V1(1,:),'r')
title('Vin vs Vout with Noise (Standard Cap)')
xlabel( 'Time (ms)')
ylabel('Voltage')

%Fourier
freq = linspace(-500,500,steps);

FV1 = fft(V1.');

FSV1 = fftshift(FV1);

figure(2)
plot(freq, 20*log10(abs(FSV1(:,1))),'r')
hold on 
plot(freq, 20*log10(abs(FSV1(:,7))),'b')
title('Fourier Transform (Standard Cap)')
axis([-500 500 -100 inf])
ylabel('Voltage')
xlabel('Frequency')


%% Changing Capcitance Values

% Using C2 (Large Cap)

V1 = zeros(7,1);
Gpulse = zeros(7,1);
Vprev = zeros(7,1);


for i = 1:steps
    
    % Current Source
    In = 0.001*randn();
    
    
    Gpulse(4,1) = In;
    
    Gpulse(1,1) = exp(-1/2*((i/steps-0.06)/(0.03))^2);
         
    V1(:,i) = (C2./dt+G)\(Gpulse+C2*Vprev/dt);
    
    Vprev = V1(:,i);
    
end

figure(3)
plot(Xval, V1(7,:),'b')
hold on
plot(Xval, V1(1,:),'r')
title('Vin vs Vout with Noise (Large Cap)')
xlabel( 'Time (ms)')
ylabel('Voltage')

%Fourier
freq = linspace(-500,500,steps);

FV1 = fft(V1.');

FSV1 = fftshift(FV1);

figure(4)
plot(freq, 20*log10(abs(FSV1(:,1))),'r')
hold on 
plot(freq, 20*log10(abs(FSV1(:,7))),'b')
axis([-500 500 -100 inf])
title('Fourier Transform (Large Cap)')
ylabel('Voltage')
xlabel('Frequency')

%% Using C3 (Small Cap)

V1 = zeros(7,1);
Gpulse = zeros(7,1);
Vprev = zeros(7,1);


for i = 1:steps
    
    % Current Source
    In = 0.001*randn();
    
    
    Gpulse(4,1) = In;
    
    Gpulse(1,1) = exp(-1/2*((i/steps-0.06)/(0.03))^2);
         
    V1(:,i) = (C3./dt+G)\(Gpulse+C3*Vprev/dt);
    
    Vprev = V1(:,i);
    
end




figure(5)
plot(Xval, V1(7,:),'b')
hold on
plot(Xval, V1(1,:),'r')
title('Vin vs Vout with Noise (Small Cap)')
xlabel( 'Time (ms)')
ylabel('Voltage')

%Fourier
freq = linspace(-500,500,steps);

FV1 = fft(V1.');


FSV1 = fftshift(FV1);



figure(6)
plot(freq, 20*log10(abs(FSV1(:,1))),'r')
hold on 
plot(freq, 20*log10(abs(FSV1(:,7))),'b')
title('Fourier Transform (Small Cap)')
axis([-500 500 -100 inf])
ylabel('Voltage')
xlabel('Frequency')

%% Changing Time Step:

% Large Time Step
steps2 = 500;
Xval = 1:steps2;
dt = 1/steps2;

V1 = zeros(7,1);
Gpulse = zeros(7,1);
Vprev = zeros(7,1);




for i = 1:steps2
    
    % Current Source
    In = 0.001*randn();
    
    
    Gpulse(4,1) = In;
    
    Gpulse(1,1) = exp(-1/2*((i/steps2-0.06)/(0.03))^2);
         
    V1(:,i) = (C1./dt+G)\(Gpulse+C1*Vprev/dt);
    
    Vprev = V1(:,i);
    
end

figure(7)
plot(Xval, V1(7,:),'b')
hold on
plot(Xval, V1(1,:),'r')
title('Vin vs Vout with Noise (Large Tstep)')
xlabel( 'Time (ms)')
ylabel('Voltage')

%Fourier
freq = linspace(-500,500,steps2);

FV1 = fft(V1.');

FSV1 = fftshift(FV1);

figure(8)
plot(freq, 20*log10(abs(FSV1(:,1))),'r')
hold on 
plot(freq, 20*log10(abs(FSV1(:,7))),'b')
title('Fourier Transform (Large Tstep)')
axis([-500 500 -100 inf])
xlabel( 'Frequency')
ylabel('Voltage')

%% Small Time Step 

steps3 = 2500;
Xval = 1:steps3;
dt = 1/steps3;

V1 = zeros(7,1);
Gpulse = zeros(7,1);
Vprev = zeros(7,1);


for i = 1:steps3
    
    % Current Source
    In = 0.001*randn();
    
    
    Gpulse(4,1) = In;
    
    Gpulse(1,1) = exp(-1/2*((i/steps3-0.06)/(0.03))^2);
         
    V1(:,i) = (C1./dt+G)\(Gpulse+C1*Vprev/dt);
    
    Vprev = V1(:,i);
    
end

figure(9)
plot(Xval, V1(7,:),'b')
hold on
plot(Xval, V1(1,:),'r')
title('Vin vs Vout with Noise (Small Tstep)')
xlabel( 'Time (ms)')
ylabel('Voltage')

%Fourier
freq = linspace(-500,500,steps3);

FV1 = fft(V1.');

FSV1 = fftshift(FV1);

figure(10)
plot(freq, 20*log10(abs(FSV1(:,1))),'r')
hold on 
plot(freq, 20*log10(abs(FSV1(:,7))),'b')
title('Fourier Transform (Small Tstep)')
axis([-500 500 -100 inf])
xlabel( 'Frequency')
ylabel('Voltage')