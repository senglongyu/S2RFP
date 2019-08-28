clear;
clc;
close all;

format long;

file1 = '8153-25.4mm.s2p';
file2 = '8153-10mm.s2p';

% CPW length from file1
Li = 25.4e-3;
% CPW length from file2
Lj = 10e-3;

[f1,S1,freq_noise1,data_noise1,Zo1] = SXPParse(file1);
[f2,S2,freq_noise2,data_noise2,Zo2] = SXPParse(file2);

% Check if input vectors are same length and same content
if (length(f1) ~= length(f2))
    error("Input vectors are not the same size! Hint: check SnP files vector lengths.");
elseif (f1 ~= f2)
    error("Frequency contents in SnP files don't match!");
end

f = f1;
n = length(f);

% Initialization
Mi(1:2,1:2,1:n) = zeros;
Mj(1:2,1:2,1:n) = zeros;
Mij(1:2,1:2,1:n) = zeros;
lambdaij_1M(1:2,1:2,1:n) = zeros;
lambdaij_2M(1:2,1:2,1:n) = zeros;
lambdaij(1:2,1:2,1:n) = zeros;
gamma(1:n) = zeros;
alpha(1:n) = zeros;
beta(1:n) = zeros;

% Matrix operations
for x = 1:n
    Mi(1,1,x) = 1/S1(2,1,x)*(S1(1,2,x)*S1(2,1,x)-S1(1,1,x)*S1(2,2,x));
    Mi(1,2,x) = 1/S1(2,1,x)*S1(1,1,x);
    Mi(2,1,x) = 1/S1(2,1,x)*(-S1(2,2,x));
    Mi(2,2,x) = 1/S1(2,1,x);
    
    Mj(1,1,x) = 1/S2(2,1,x)*(S2(1,2,x)*S2(2,1,x)-S2(1,1,x)*S2(2,2,x));
    Mj(1,2,x) = 1/S2(2,1,x)*S2(1,1,x);
    Mj(2,1,x) = 1/S2(2,1,x)*(-S2(2,2,x));
    Mj(2,2,x) = 1/S2(2,1,x);
    
    Mij(:,:,x) = Mj(:,:,x)*inv((Mi(:,:,x)));
    
    lambdaij_1M(x) = ((Mij(1,1,x)+Mij(2,2,x))+sqrt((Mij(1,1,x)-Mij(2,2,x))^2+4*Mij(1,2,x)*Mij(2,1,x)))/2;
    lambdaij_2M(x) = ((Mij(1,1,x)+Mij(2,2,x))-sqrt((Mij(1,1,x)-Mij(2,2,x))^2+4*Mij(1,2,x)*Mij(2,1,x)))/2;
    
    lambdaij(x) = (1/2)*(lambdaij_1M(x)+(1/lambdaij_2M(x)));
    
    gamma(x) = log(lambdaij(x))/(Li-Lj);
    
    % Force alpha to be positive
    if (real(gamma(x))<0)
        gamma(x) = -gamma(x);
    end
    
    alpha(x) = real(gamma(x));
    beta(x) = imag(gamma(x));
    
end

% Phase Correction
% beta = beta*(Li-Lj);
% for x = 2:n
%     tol = 0.5; % Tolerance in radians
%     delta = diff([beta(x-1),beta(x)]);
%     if delta > tol
%         beta(x-1) = beta(x-1) - delta;
%     end
% end
% beta = beta/(Li-Lj);

% Unwrap beta
% First scale back to [-pi, pi] range and scale back up
beta_w = beta;
beta_s = beta*(Li-Lj);
beta_s = unwrap(beta_s);
beta = beta_s/(Li-Lj);
gamma(:) = alpha(:) +1i*beta(:);

h1 = 1.54e-3;  % lower substrate height
h2 = 71e-6; % sample height
er1 = 2.33; % lower substrate er
w = 1e-3;
s = 0.5e-3;
ur = 1;
u0 = 4*pi*1e-7;
u = ur*u0;
e0 = 8.854e-12;

k0 = w/(w+2*s);
k1 = sinh(pi*w/(4*h1))/sinh((pi*(w+2*s))/(4*h1));
k2 = sinh(pi*w/(4*h2))/sinh((pi*(w+2*s))/(4*h2));
k0p = sqrt(1-k0^2);
k1p = sqrt(1-k1^2);
k2p = sqrt(1-k2^2);

q1 = 0.5*(ellipke(k1)/ellipke(k1p))/(ellipke(k0)/ellipke(k0p));
q2 = 0.5*(ellipke(k2)/ellipke(k2p))/(ellipke(k0)/ellipke(k0p));

% Initialization
e(1:n) = zeros;
eeff(1:n) = zeros;
er(1:n) = zeros;
tand(1:n) = zeros;

for x = 1:n
    e(x) = (alpha(x)^2-beta(x)^2)/(-(2*pi*f(x))^2*u);
    eeff(x) = e(x)/e0;
    er(x) = (eeff(x)-q1*(er1-1)+q2*er1-1)/q2; % conversion to er, see notes
    tand(x) = 2*alpha(x)*beta(x)/((2*pi*f(x))^2*u*e(x));
end

for x = 1:n
    k = (((S1(1,1,x)^2-S1(2,1,x)^2+1)^2-(2*S1(1,1,x))^2)/((2*S1(2,1,x))^2))^(1/2);
    EE = 1/((1-S1(1,1,x)^2+S1(2,1,x)^2)/(2*S1(2,1,x))+k);
    gamma1(x) = -(1/Li)*log(EE);
    if real(gamma1(x))<0
        EE = 1/((1-S1(1,1,x)^2+S1(2,1,x)^2)/(2*S1(2,1,x))-k);
        gamma1(x) = -(1/Li)*log(EE);
    end
    
    alpha1(x) = real(gamma1(x));
    beta1(x) = imag(gamma1(x));
end

beta1_w = beta1;
beta1_s = beta1*Li;
beta1_s = unwrap(beta1_s);
beta1 = beta1_s/Li;
gamma1(:) = alpha1(:)+1i*beta1(:);

for x = 1:n
    es(x) = (alpha1(x)^2-beta1(x)^2)/(-(2*pi*f(x))^2*u);
    eeffs(x) = es(x)/e0;
    ers(x) = (eeffs(x)-q1*(er1-1)+q2*er1-1)/q2; % conversion to er, see notes
    tands(x) = 2*alpha1(x)*beta1(x)/((2*pi*f(x))^2*u*es(x));
end
f_start = 1; % 159;
f_end = length(f); %1631;

figure(1)
plot(f/1e9,abs(er),f/1e9,ers,'--');
xlim([f(f_start)/1e9 f(f_end)/1e9]);
ylim([25 35]);
xlabel('Frequency (GHz)');
ylabel("Dielectric Constant (\epsilon_r)");
legend('Method 1','Method 2')

figure(2)
plot(f/1e9,tand,f/1e9,tands,'--');
xlim([f(f_start)/1e9 f(f_end)/1e9]);
ylim([0 0.1]);
xlabel('Frequency (GHz)');
ylabel("Loss Tangent (tan\delta)");
legend('Method 1','Method 2')
figure(3);

subplot(2,1,1);
plot(f/1e9,alpha,f/1e9,alpha1,'--');
xlim([f(f_start)/1e9 f(f_end)/1e9]);
xlabel('Frequency (GHz)');
title('Attenuation Constant \alpha');
legend('Method 1','Method 2')
subplot(2,1,2);
plot(f/1e9,beta,f/1e9,beta1,'--');
xlim([f(f_start)/1e9 f(f_end)/1e9]);
xlabel('Frequency (GHz)');
title('Phase Constant \beta');
legend('Method 1','Method 2')

figure(4);
eri = 5;
error1 = abs(er-eri);
error2 = abs(ers-eri);
plot(f/1e9,error1,f/1e9,error2,'--');
ylim([0 10]);
xlabel('Frequency (GHz)');
ylabel("Dielectric Constant Error");
legend('Method 1','Method 2')

figure(5);
ELi = wrapTo180(beta*Li*180/pi);
ELj = wrapTo180(beta*Lj*180/pi);
plot(f/1e9,ELi,f/1e9,ELj,'--');
ylabel('Electrical Length (°)');
xlabel('Frequency (GHz)');
legend('25.4 mm','10 mm');


figure(6);
for x = 1:n
    s11_1(x) = S1(1,1,x);
    s11_2(x) = S2(1,1,x);
    s21_1(x) = S1(2,1,x);
    s21_2(x) = S2(2,1,x);
end
sp1_angle = angle(s21_1);
sp2_angle = angle(s21_2);
sp1 = abs(s11_1);
sp2 = abs(s11_2);
subplot(2,1,1);
plot(f/1e9,20*log10(sp1),f/1e9,20*log10(sp2),'--');
subplot(2,1,2);
plot(f/1e9,sp1_angle,f/1e9,sp2_angle,'--');

% 
% disp(mean(er(f_start:f_end)));
% disp(mean(tand(f_start:f_end)));
% disp(mean(ers(f_start:f_end)));
% disp(mean(tands(f_start:f_end)));
% figure(1)
% subplot(2,2,1);
% plot(f,er);
% xlim([f(1) f(n)]);
% xlabel('Frequency (Hz)');
% title("Dielectric Constant");
% subplot(2,2,2);
% plot(f,tand);
% xlim([f(1) f(n)]);
% xlabel('Frequency (Hz)');
% title("Loss Tangent");
% subplot(2,2,3);
% plot(f,alpha);
% xlim([f(1) f(n)]);
% xlabel('Frequency (Hz)');
% title("Alpha");
% subplot(2,2,4);
% plot(f,beta);
% xlim([f(1) f(n)]);
% xlabel('Frequency (Hz)');
% title("Beta");