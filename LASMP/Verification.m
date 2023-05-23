clc;
clear;

%
% Define global parameters
%
Mu = 1.0e6;  %Shear modulus of initial phase
Mu1 = 3.0e6; %Shear modulus of light-induced phase
ef = 3;      %Stretch at which light is induced
af = 0.8;    %Final extent of reaction
kf = 0.0025; %Reaction constant of crosslinking
kr = 0.0085; %Reaction constant of clevage

tf = 1000;   %Time of crosslinking
tu = 1000;   %Time of unloading
tr = 1000;   %Time of clevage

%
% True strain vs. True stress
%
figure(1)
hold on;
%%--Loading--%%
lamda = 1:0.1:ef;
T1 = lamda.^2-1./lamda;

plot((lamda-1)*100,T1,'k');


%%--Unloading--%%
lamda = ef:-0.1:1;
T1 = lamda.^2-1./lamda;
T2 = 0.7145*(Mu1/Mu)*(lamda.^2/ef^2-ef./lamda);  
T = T1 + T2;

plot((lamda-1)*100,T1,'k',(lamda-1)*100,T,'k','linewidth',1.2);

%%--Simulation Results--%%
[strain stress] = textread('NominalStrainVsTrueStress.txt','%f%f','headerlines',3);
plot(strain(1:3:end)*100,stress(1:3:end)/Mu,'ko');

axis([0 250 0 12]);
box on;
xlabel('Engineering Strain (%)');
ylabel('Normalized Cauchy Stress');
legend('Analytical solution','','Numerical solution');
legend boxoff;
set(gca,'FontSize',12);
hold off;



%
% Time vs. Extent of reaction
%
figure(2);
hold on;

dt = 1;
alpha = zeros(1,max(size(0:dt:4000)));
alpha(1) = 0;
time = zeros(1,max(size(0:dt:4000)));
time(1) = 0;

for i = 2:max(size(0:dt:4000));
    time(i) = time(i-1) + dt;  
    if (time(i) >0 && time(i) <= 1000)
        alpha(i) = 0;
    elseif (time(i) > 1000 && time(i) <= 2000)
        alpha(i) = alpha(i-1) + dt*kf*(1-alpha(i-1))^2;
        if(alpha(i) >= af)
            alpha(i) = af;
        end
    elseif (time(i) > 2000 && time(i) <= 3000)
        alpha(i) = alpha(i-1);
    else
        alpha(i) = alpha(i-1) - dt*kr*alpha(i-1);
    end
end 
plot(time,alpha,'k','linewidth',1.2);

%--Simulation Results--%%
[t a] = textread('ExtentofReaction.txt','%f%f','headerlines',3);
plot(t(1:3:end),a(1:3:end),'ko');

axis([0 4000 0 1]);
box on;
xlabel('Time (s)');
ylabel('Extent of Reaction');
legend('Analytical solution','Numerical solution');
legend boxoff;
set(gca,'FontSize',12);
hold off;
