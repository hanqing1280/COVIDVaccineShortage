%----------------------------Introduction------------------------------%
% South Africa Data Fit
% Method: MCMC WITHOUT Initial Values Being Parameters
% Data Set: Daily New Cases + Daily New Deaths
% Input: None
%----------------------------------------------------------------------%

clear model data params options
load('AfricaCovidData');

global np nps sigma gammai gammaa c nsimu paramstr nparam color;
np=58558000*[0.3661;0.4095;0.1635;0.0609];
nps=sum(np);npperc=np./nps;
sigma=1./5.1;gammai=1./4.2;gammaa=1./5.88;
c=[10 8 8 8;
   8 10 8 8;
   8 8 8 8;
   8 8 8 8];
nsimu=1000;
paramstr={'\beta^a_1';'\beta^a_2';'\beta^a_3';'\beta^a_4';'\rho^a_i';'p^a_i';'\delta^a_1';'\delta^a_2';'\delta^a_3';'\delta^a_4';'\eta^a_1';'\eta^a_2';'\eta^a_3';'\eta^a_4';'\phi^a_1';'\phi^a_2';'\phi^a_3';'\phi^a_4';'q^a_1';'q^a_2';'q^a_3';'q^a_4'};
nparam=length(paramstr);

%colors
color=[252 41 30;250 200 205;219 249 244;54 195 201;0 70 222;255 255 199;255 170 50]./255;

data.ydata=[SouthAfricaNewCase(250:368) SouthAfricaNewDeath(228:346)];
data.time=(0:length(data.ydata(:,1))-1)';
data.ydata=[data.time data.ydata];
e0=[552.939,612.716,239.109,89.6658]'*100./sigma;
d0=[27.1173,30.0489,11.7264,4.3974]';
s0=np-e0-d0;
data.y0=[s0;e0;ones(4,1);ones(4,1);ones(4,1);d0;ones(4,1)];

% calling fminsearch
theta00=[0.5,0.5,0.5,0.5,0.55,0.792,0.244,0.244,0.244,0.244,0.5,0.5,0.5,0.5,0.1,0.1,0.1,0.1,0.032,0.032,0.032,0.032]';
[theta0,ss0]=fminsearch(@strainAmodelss,theta00,[],data)
mse = ss0/(length(data.ydata)-nparam);

params = {
   %  name,      init,       min,   max
    {'\beta_1',    theta0(1),   0,     1}
    {'\beta_2',    theta0(2),   0,     1}
    {'\beta_3',    theta0(3),   0,     1}
    {'\beta_4',    theta0(4),   0,     1}
    {'\rho_i',      theta0(5),  0,     1}
    {'p_i' ,       theta0(6),   0,     1}
    {'\delta_1',   theta0(7),   0,     1}
    {'\delta_2',   theta0(8),   0,     1}
    {'\delta_3',   theta0(9),   0,     1}
    {'\delta_4',   theta0(10),  0,     1}
    {'\eta_1',     theta0(11),  0,     1}
    {'\eta_2',     theta0(12),  0,     1}
    {'\eta_3',     theta0(13),  0,     1}
    {'\eta_4',     theta0(14),  0,     1}
    {'\phi_1',     theta0(15),  0,     1}
    {'\phi_2',     theta0(16),  0,     1}
    {'\phi_3',     theta0(17),  0,     1}
    {'\phi_4',     theta0(18),  0,     1}
    {'q_1',       theta0(19),   0,     1}
    {'q_2',       theta0(20),   0,     1}
    {'q_3',       theta0(21),   0,     1}
    {'q_4',       theta0(22),   0,     1}
         };
         
model.ssfun = @strainAmodelss;
model.sigma2 = mse;

options.nsimu = nsimu;
options.updatesigma = 1; 
options.method= 'dram'; % adaptation method, 'mh', 'dr', 'am', or 'dram'
options.adaptint= 500;    % how often to adapt the proposal

[results,chain,s2chain] = mcmcrun(model,data,params,options);

chainstats(chain,results)

%-------------------------------Figures------------------------------------%
%-----------Parameter Chain Plot-----------%
clf
close all
figure(1)
mcmcplot(chain,[],results,'chainpanel')
%mcmcplot(sqrt(s2chain),[],[],'dens',2)
%title('Error Std')
sgtitle('South Africa');

%-------Parameter Histogram Plot--------%
figure(2);
for i=1:nparam
    subplot(5,8,i)
    histogram(chain(:,i),nsimu);
    title(paramstr{i},'FontSize',20);
end
sgtitle('South Africa');

%---------------------Fit Plot with MCMC Mean---------------------%
figure(3)
[t,cases,deaths]=strainAcasedeath(data.time,mean(chain),data.y0);
plot(data.time,data.ydata(:,2),'o','Color',color(1,:),'LineWidth',2)
hold on
plot(t,cases,'-','Color',color(1,:),'LineWidth',2)
hold on
plot(data.time,data.ydata(:,3),'o','Color',color(4,:),'LineWidth',2)
hold on
plot(t,deaths,'-','Color',color(4,:),'LineWidth',2)
legend('Data of new cases','Model of new cases','Data of deaths','Model of deaths','Location','north west');
title('MCMC Fit (Fixed Init. Values)','FontSize',42);
xlabel('Day','FontSize',30);ylabel('Number','FontSize',30);

%------------------Fit Plot with LSQ Estimators-------------------%
figure(4)
[t,cases,deaths]=strainAcasedeath(data.time,theta0,data.y0);
plot(data.time,data.ydata(:,2),'o','Color',color(1,:),'LineWidth',2)
hold on
plot(t,cases,'-','Color',color(1,:),'LineWidth',2)
hold on
plot(data.time,data.ydata(:,3),'o','Color',color(4,:),'LineWidth',2)
hold on
plot(t,deaths,'-','Color',color(4,:),'LineWidth',2)
legend('Data of new cases','Model of new cases','Data of new deaths','Model of new deaths','Location','north west');
xlabel('Day','FontSize',30);ylabel('Number','FontSize',30);
title('LSQ Fit (Fixed Init. Values)','FontSize',42);

%------------------Predictive Envelop Plot-------------------%
%figure(5)
%out = mcmcpred(results,chain,[],data,@strainAcase);
%mcmcpredplot(out,2,1);
%hold on
%plot(data.time,data.ydata(:,2),'s','Color',color(2,:),'LineWidth',2); 
%xlabel('Day','FontSize',30);ylabel('Number','FontSize',30);
%title('MCMC Fit','FontSize',42);
%set(gcf,'unit','centimeters','position',[10 15 20 15])

%-------------------------------Functions------------------------------------%
function casemod=strainAcase(theta,data)
global sigma;
[t,y]=strainAmodel(data.time,theta,data.y0);
casemod=theta(6)*sigma*(y(:,5)+y(:,6)+y(:,7)+y(:,8));
end  
  

function ss=strainAmodelss(theta,data)
time=data.time; 
y0=data.y0;
[t,cases,deaths]=strainAcasedeath(time,theta,y0);
ss=sum((data.ydata(:,2)-cases).^2);%+sum((data.ydata(:,3)-deaths).^2); % the total SS
end


function [t,cases,deaths]=strainAcasedeath(time,theta,y0)
global sigma;
[t,y]=strainAmodel(time,theta,y0);
cases=theta(6)*sigma*(y(:,5)+y(:,6)+y(:,7)+y(:,8));
deaths=theta(7)*y(:,9)+theta(8)*y(:,10)+theta(9)*y(:,11)+theta(10)*y(:,12)+theta(19)*theta(15)*y(:,17)+theta(20)*theta(16)*y(:,18)+theta(21)*theta(17)*y(:,19)+theta(22)*theta(18)*y(:,20);
end  


function [t,y]=strainAmodel(time,theta,y0)
[t,y]=ode45(@strainAode,time,y0,[],theta);
end


function dy=strainAode(t,y,theta)
% known parameter values
global np sigma gammai gammaa c;

% take parameters and components out from y and theta
beta1=theta(1);beta2=theta(2);beta3=theta(3);beta4=theta(4);
rho=theta(5);p=theta(6);
delta1=theta(7);delta2=theta(8);delta3=theta(9);delta4=theta(10);
eta1=theta(11);eta2=theta(12);eta3=theta(13);eta4=theta(14);
phi1=theta(15);phi2=theta(16);phi3=theta(17);phi4=theta(18);
q1=theta(19);q2=theta(20);q3=theta(21);q4=theta(22);

% variables
s1=y(1);s2=y(2);s3=y(3);s4=y(4);
e1=y(5);e2=y(6);e3=y(7);e4=y(8);
i1=y(9);i2=y(10);i3=y(11);i4=y(12);
a1=y(13);a2=y(14);a3=y(15);a4=y(16);
h1=y(17);h2=y(18);h3=y(19);h4=y(20);
d1=y(21);d2=y(22);d3=y(23);d4=y(24);
r1=y(25);r2=y(26);r3=y(27);r4=y(28);

% define the ODE
dy(1)=-beta1*s1*(c(1,1)*(i1+rho*a1)/np(1)+c(1,2)*(i2+rho*a2)/np(2)+c(1,3)*(i3+rho*a3)/np(3)+c(1,4)*(i4+rho*a4)/np(4));
dy(2)=-beta2*s2*(c(2,1)*(i1+rho*a1)/np(1)+c(2,2)*(i2+rho*a2)/np(2)+c(2,3)*(i3+rho*a3)/np(3)+c(2,4)*(i4+rho*a4)/np(4));
dy(3)=-beta3*s3*(c(3,1)*(i1+rho*a1)/np(1)+c(3,2)*(i2+rho*a2)/np(2)+c(3,3)*(i3+rho*a3)/np(3)+c(3,4)*(i4+rho*a4)/np(4));
dy(4)=-beta4*s4*(c(4,1)*(i1+rho*a1)/np(1)+c(4,2)*(i2+rho*a2)/np(2)+c(4,3)*(i3+rho*a3)/np(3)+c(4,4)*(i4+rho*a4)/np(4));
dy(5)=beta1*s1*(c(1,1)*(i1+rho*a1)/np(1)+c(1,2)*(i2+rho*a2)/np(2)+c(1,3)*(i3+rho*a3)/np(3)+c(1,4)*(i4+rho*a4)/np(4))-sigma*e1;
dy(6)=beta2*s2*(c(2,1)*(i1+rho*a1)/np(1)+c(2,2)*(i2+rho*a2)/np(2)+c(2,3)*(i3+rho*a3)/np(3)+c(2,4)*(i4+rho*a4)/np(4))-sigma*e2;
dy(7)=beta3*s3*(c(3,1)*(i1+rho*a1)/np(1)+c(3,2)*(i2+rho*a2)/np(2)+c(3,3)*(i3+rho*a3)/np(3)+c(3,4)*(i4+rho*a4)/np(4))-sigma*e3;
dy(8)=beta4*s4*(c(4,1)*(i1+rho*a1)/np(1)+c(4,2)*(i2+rho*a2)/np(2)+c(4,3)*(i3+rho*a3)/np(3)+c(4,4)*(i4+rho*a4)/np(4))-sigma*e4;
dy(9)=p*sigma*e1-(delta1+eta1+gammai)*i1;
dy(10)=p*sigma*e2-(delta2+eta2+gammai)*i2;
dy(11)=p*sigma*e3-(delta3+eta3+gammai)*i3;
dy(12)=p*sigma*e4-(delta4+eta4+gammai)*i4;
dy(13)=(1-p)*sigma*e1-gammaa*a1;
dy(14)=(1-p)*sigma*e2-gammaa*a2;
dy(15)=(1-p)*sigma*e3-gammaa*a3;
dy(16)=(1-p)*sigma*e4-gammaa*a4;
dy(17)=eta1*i1-phi1*h1;
dy(18)=eta2*i2-phi2*h2;
dy(19)=eta3*i3-phi3*h3;
dy(20)=eta4*i4-phi4*h4;
dy(21)=delta1*i1+q1*phi1*h1;
dy(22)=delta2*i2+q2*phi2*h2;
dy(23)=delta3*i3+q3*phi3*h3;
dy(24)=delta4*i4+q4*phi4*h4;
dy(25)=gammaa*a1+gammai*i1+(1-q1)*phi1*h1;
dy(26)=gammaa*a2+gammai*i2+(1-q2)*phi2*h2;
dy(27)=gammaa*a3+gammai*i3+(1-q3)*phi3*h3;
dy(28)=gammaa*a4+gammai*i4+(1-q4)*phi4*h4;

dy=dy(:);% make sure that we return a column vector
end


