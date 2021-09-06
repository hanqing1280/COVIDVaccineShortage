%----------------------------Introduction------------------------------%
% South Africa Data Fit
% Method: MCMC WITH Initial Values Being Parameters
% Data Set: Daily New Cases
% Input: None
%----------------------------------------------------------------------%
addpath('/Users/Qing/Desktop/IDRC/Vaccine shortage/Vaccine shortage MATLAB Codes/mcmcstat-master')

clear model data params options
load('AfricaCovidData');

global np nps npperc sigma gammai gammaa c nsimu paramstr nparam time color country h0 d0 e0;
country='SouthAfrica';
np=59308690*[0.3661;0.4095;0.1635;0.0609];%population structure
nps=sum(np);npperc=np./nps;
sigma=1./5.1;gammai=1./4.2;gammaa=1./5.88;
c=[10 8 8 8;
   8 10 8 8;
   8 8 8 8;
   8 8 8 8];
nsimu=1000;
paramstr={'\beta^a_1';'\beta^a_2';'\beta^a_3';'\beta^a_4';'\rho^a_i';'p^a_i';'\delta^a_1';'\delta^a_2';'\delta^a_3';'\delta^a_4';'\eta^a_1';'\eta^a_2';'\eta^a_3';'\eta^a_4';'\phi^a_1';'\phi^a_2';'\phi^a_3';'\phi^a_4';'q^a_1';'q^a_2';'q^a_3';'q^a_4';'A^a_{10}';'A^a_{20}';'A^a_{30}';'A^a_{40}';'R^a_{10}';'R^a_{20}';'R^a_{30}';'R^a_{40}'};
nparam=length(paramstr);

%colors
color=[252 41 30;250 200 205;219 249 244;54 195 201;0 70 222;255 255 199;255 170 50]./255;

coviddata=[SouthAfricaNewCase(250:368) SouthAfricaNewDeath(228:346)];
time=(0:length(coviddata(:,1))-1)';
data.ydata=[time coviddata];
h0=4848;
d0=SouthAfricaTotalDeath(228);
e0=SouthAfricaNewCase(250);


% calling fminsearch
theta00=[0.5,0.5,0.5,0.5,0.55,0.792,0.0244,0.0244,0.0244,0.0244,0.5,0.5,0.5,0.5,0.1,0.1,0.1,0.1,0.032,0.032,0.032,0.032,100*ones(1,8)]';
[theta0,ss0]=fminsearch(@strainAmodelss,theta00,[],data)
mse = ss0/(length(data.ydata)-nparam);

params = {
   %  name,      init,       min,   max
    {'\beta^a_1',    theta00(1),   0,     1}
    {'\beta^a_2',    theta00(2),   0,     1}
    {'\beta^a_3',    theta00(3),   0,     1}
    {'\beta^a_4',    theta00(4),   0,     1}
    {'\rho^a_i',     theta00(5),   0,     1}
    {'p^a_i' ,       theta00(6),   0,     1}
    {'\delta^a_1',   theta00(7),   0,     1}
    {'\delta^a_2',   theta00(8),   0,     1}
    {'\delta^a_3',   theta00(9),   0,     1}
    {'\delta^a_4',   theta00(10),  0,     1}
    {'\eta^a_1',     theta00(11),  0,     1}
    {'\eta^a_2',     theta00(12),  0,     1}
    {'\eta^a_3',     theta00(13),  0,     1}
    {'\eta^a_4',     theta00(14),  0,     1}
    {'\phi^a_1',     theta00(15),  0,     1}
    {'\phi^a_2',     theta00(16),  0,     1}
    {'\phi^a_3',     theta00(17),  0,     1}
    {'\phi^a_4',     theta00(18),  0,     1}
    {'q^a_1',        theta00(19),  0,     1}
    {'q^a_2',        theta00(20),  0,     1}
    {'q^a_3',        theta00(21),  0,     1}
    {'q^a_4',        theta00(22),  0,     1}
    {'A^a_{10}',     theta00(23),  0,     nps}
    {'A^a_{20}',     theta00(24),  0,     nps}
    {'A^a_{30}',     theta00(25),  0,     nps}
    {'A^a_{40}',     theta00(26),  0,     nps}
    {'R^a_{10}',     theta00(27),  0,     nps}
    {'R^a_{20}',     theta00(28),  0,     nps}
    {'R^a_{30}',     theta00(29),  0,     nps}
    {'R^a_{40}',     theta00(30),  0,     nps}
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
sgtitle('South Africa','FontSize',30);
saveas(gcf,['/Users/Qing/Desktop/IDRC/Vaccine shortage/Latex/',[country,'ParaChainOnlyNewCase.jpg']]);

%-------Parameter Histogram Plot--------%
figure(2);
for i=1:nparam
    subplot(3,10,i)
    histogram(chain(:,i),nsimu);
    title(paramstr{i},'FontSize',20);
end
sgtitle('South Africa','FontSize',30);
saveas(gcf,['/Users/Qing/Desktop/IDRC/Vaccine shortage/Latex/',[country,'ParaDensityOnlyNewCase.jpg']]);

%---------------------Fit Plot with MCMC Mean---------------------%
figure(3)
[t,cases,deaths]=strainAcasedeath(mean(chain));
plot(time,data.ydata(:,2),'o','Color',color(1,:),'LineWidth',2)
hold on
plot(t,cases,'-','Color',color(1,:),'LineWidth',2)
hold on
plot(time,data.ydata(:,3),'o','Color',color(4,:),'LineWidth',2)
hold on
plot(t,deaths,'-','Color',color(4,:),'LineWidth',2)
legend('Data of new cases','Model of new cases','Data of deaths','Model of deaths','Location','north west');
title('MCMC Fit','FontSize',42);
xlabel('Day','FontSize',30);ylabel('Number','FontSize',30);
saveas(gcf,['/Users/Qing/Desktop/IDRC/Vaccine shortage/Latex/',[country,'MCMCFitOnlyNewCase.jpg']]);


%------------------Fit Plot with LSQ Estimators-------------------%
figure(4)
[t,cases,deaths]=strainAcasedeath(theta0);
plot(time,data.ydata(:,2),'o','Color',color(1,:),'LineWidth',2)
hold on
plot(t,cases,'-','Color',color(1,:),'LineWidth',2)
hold on
plot(time,data.ydata(:,3),'o','Color',color(4,:),'LineWidth',2)
hold on
plot(t,deaths,'-','Color',color(4,:),'LineWidth',2)
legend('Data of new cases','Model of new cases','Data of new deaths','Model of new deaths','Location','north west');
xlabel('Day','FontSize',30);ylabel('Number','FontSize',30);
title('LSQ Fit','FontSize',42);
saveas(gcf,['/Users/Qing/Desktop/IDRC/Vaccine shortage/Latex/',[country,'LSQFitOnlyNewCase.jpg']]);


%------------------Predictive Envelop Plot-------------------%
%figure(5)
%out = mcmcpred(results,chain,[],data,@strainAcase,500);
%mcmcpredplot(out,2,1);
%hold on
%plot(time,data.ydata(:,2),'s','Color',color(1,:),'LineWidth',2); 
%xlabel('Day','FontSize',30);ylabel('Number','FontSize',30);
%title('MCMC Fit','FontSize',42);

%-------------------------------Functions------------------------------------%
function casemod=strainAcase(theta)
global sigma time np;
[t,y]=strainAmodel(time,theta);
casemod=theta(6)*sigma*(y(:,5)*np(1)+y(:,6)*np(2)+y(:,7)*np(3)+y(:,8)*np(4));
end  
   

function ss=strainAmodelss(theta,data)
[t,cases,deaths]=strainAcasedeath(theta);
ss=sum((data.ydata(:,2)-cases).^2); % the total SS
end


function [t,cases,deaths]=strainAcasedeath(theta)
global sigma time np;
[t,y]=strainAmodel(time,theta);
cases=theta(6)*sigma*(y(:,5)*np(1)+y(:,6)*np(2)+y(:,7)*np(3)+y(:,8)*np(4));
deaths=theta(7)*y(:,9)*np(1)+theta(8)*y(:,10)*np(2)+theta(9)*y(:,11)*np(3)+theta(10)*y(:,12)*np(4)+theta(19)*theta(15)*y(:,17)*np(1)+theta(20)*theta(16)*y(:,18)*np(2)+theta(21)*theta(17)*y(:,19)*np(3)+theta(22)*theta(18)*y(:,20)*np(4);
end 


function [t,y]=strainAmodel(time,theta)
global np sigma npperc h0 d0 e0;
y0=zeros(28,1);
y0(13)=theta(23)/np(1);y0(14)=theta(23)/np(2);y0(15)=theta(23)/np(3);y0(16)=theta(23)/np(4);
y0(25)=theta(27)/np(1);y0(26)=theta(28)/np(2);y0(27)=theta(29)/np(3);y0(28)=theta(30)/np(4);
%y0(13:16)=theta(23:26)./np;y0(25:28)=theta(27:30)./np;% a0;r0
y0(17)=h0*npperc(1)./np(1);y0(18)=h0*npperc(2)./np(2);
y0(19)=h0*npperc(3)./np(3);y0(20)=h0*npperc(4)./np(4);
y0(21)=d0*npperc(1)./np(1);y0(22)=d0*npperc(2)./np(2);
y0(23)=d0*npperc(3)./np(3);y0(24)=d0*npperc(4)./np(4);
%y0(17:20)=h0.*npperc./np;% h0
%y0(21:24)=d0.*npperc./np;% d0
y0(5)=e0*npperc(1)./(sigma*theta(6))./np(1);
y0(6)=e0*npperc(2)./(sigma*theta(6))./np(2);
y0(7)=e0*npperc(3)./(sigma*theta(6))./np(3);
y0(8)=e0*npperc(4)./(sigma*theta(6))./np(4);
%y0(5:8)=e0.*npperc/(sigma*theta(6))./np;% e0
y0(9)=(d0*npperc(1)-theta(19)*theta(15)*y0(17))./theta(7)./np(1);
y0(10)=(d0*npperc(2)-theta(20)*theta(16)*y0(18))./theta(8)./np(2);
y0(11)=(d0*npperc(3)-theta(21)*theta(17)*y0(19))./theta(9)./np(3);
y0(12)=(d0*npperc(4)-theta(22)*theta(18)*y0(20))./theta(10)./np(4);
%y0(9:12)=(d0.*npperc-theta(19:22).*theta(15:18).*y0(17:20))./theta(7:10)./np;% i0
y0(1)=(np(1)-(y0(5)+y0(9)+y0(13)+y0(17)+y0(21)+y0(25)))/np(1);
y0(2)=(np(2)-(y0(6)+y0(10)+y0(14)+y0(18)+y0(22)+y0(26)))/np(2);
y0(3)=(np(3)-(y0(7)+y0(11)+y0(15)+y0(19)+y0(23)+y0(27)))/np(3);
y0(4)=(np(4)-(y0(8)+y0(12)+y0(16)+y0(20)+y0(24)+y0(28)))/np(4);
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
