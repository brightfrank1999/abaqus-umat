%% Incompressible Arruda-Boyce Model 
clc;
clear;

%%Stress-Strain%%
%% True strain = ln(1 + engineering strain)
%% Ture stress = (engineering stress) * exp(true strain) = (engineering stress) * (1 + engineering strain)
N=0:0.1:3;
Lamda=zeros(max(size(N)-1),1);
EngineeringStress=zeros(max(size(N)-1),1);
TrueStress=zeros(max(size(N)-1),1);
Lamda_Bar=zeros(max(size(N))-1,1);
I_1=zeros(max(size(N))-1,1);
x=zeros(max(size(N))-1,1);

Lamda(1)=1;
EngineeringStress(1)=0;

for i=2:max(size(N)-1)
    Lamda(i)=1+i*0.1;
    I_1(i)=Lamda(i)^2+2.0/Lamda(i);
    %% Inverse Langevin(x)=(x)*(3-x^2)/(1-x^2) approximation
    Mu=2.078E6; 
    Lamda_L=2.8;   %Lamda_L=sqrt(N);
    Lamda_Bar(i)=sqrt(I_1(i)/3.0);
    x(i)=Lamda_Bar(i)/Lamda_L;
    EngineeringStress(i)=(1.0/3.0)*Mu*(3.0-x(i)^2)/(1-x(i)^2)*(Lamda(i)-1.0/Lamda(i)^2); 
    TrueStress(i)=       (1.0/3.0)*Mu*(3.0-x(i)^2)/(1-x(i)^2)*(Lamda(i)^2-1.0/Lamda(i));
end

[strain stress] = textread('EngineeringStrainVsTrueStress.txt','%f%f','headerlines',3);

%figure(1);
subplot(1,2,1)
grid on;
box on;
hold on;
plot((Lamda-1)*100,EngineeringStress,'-.','linewidth',1.4);
% plot([300,300],[0,3e7],':');
plot(strain*100,stress./(strain+1),'bo');
axis([0,380 0 2.5e7]);
xlabel('Engineering Strain (%)');
ylabel('Engineering Stress (Pa)');
legend('Analytical Results', 'Simulation Results')
legend boxoff; 

%figure(2);
subplot(1,2,2)
grid on
box on
hold on;
plot((Lamda-1)*100,TrueStress,'-.','linewidth',1.4);
plot(strain*100,stress,'bo');
axis([0,380 0 10.0e7]);
xlabel('Engineering Strain (%)');
ylabel('True Stress (Pa)');
legend('Analytical Results', 'Simulation Results')
legend boxoff; 











