clear all
%close all
clc

qsave = 0;

k      = 0.7;
phi_p  = 1;
phi_m  = 0;
Rinit  = 1;

kcrit = (phi_p+phi_m)/(2-(phi_p+phi_m))

Rss   = (3*k*(1 + k) + sqrt(3)*sqrt(k*(1 + k)^2*(4 + 3*k)))/(2*(1 + k)^(3/2))

vt = 0:0.4:50;
y0 = Rinit;

for i = 1:length(vt)-1

    % [t,y] = ode45(@(t,y)(y^2*(k-1)+4*gamma*sqrt(k+1)+...
    %                 2*(k+1)*y*gamma+y*(y*(k-1)-2*gamma*(k+1))*coth(sqrt(1+k)*y))...
    %                 /(4*sqrt(1+k)*y*(gamma)),[vt(i) vt(i+1)],y0);

    [t,y] = ode45(@(t,y)1/y-(k*(phi_m-1)+phi_m+(k*(phi_p-1)+phi_p)*coth(sqrt(1+k)*y))/(sqrt(1+k)*(phi_p-phi_m)),[vt(i) vt(i+1)],y0);
    
    
    y0 = y(end,1);
    vR = y0;
    
    figure(1)
    loglog(t(end),vR,'o')
    hold on

    % rin  = linspace(0,vR(end),500);
    % rout = linspace(vR(end),6*vR(end),500);
    % cin  = (4*k+(2*(vR(end)-k*vR(end)+2*(k+1)*gamma).*csch(sqrt(1+k)*vR(end)).*sinh(sqrt(1+k)*rin))./rin)/(4*(1+k));
    % cout = (2*k+(exp(sqrt(1+k)*(vR(end)-rout)).*(vR(end)-k*vR(end)-2*(1+k)*gamma) )./rout)./(2*(1+k));
    %rin  = linspace(0,vR(end),500);
    %rout = linspace(vR(end),6*vR(end),500);
    %cin  = (k*rin+vR(end)*csch(vR(end)*sqrt(1+k))*sinh(rin*sqrt(1+k)))./(rin+k*rin);
    %cout = (exp((vR(end)-rout)*sqrt(1+k))*(vR(end)-9*k*vR(end))+10*k*rout)./(10*(1+k)*rout);
    %figure(3)
    %plot(rin,cin)
    %hold on
    %plot(rout,cout)
    %hold off
     
    %pause(0.05)

end

% z    = linspace(0,5,5000);
% %dRdt = (z.^2.*(k-1)+4*gamma*sqrt(k+1)+...
% %                2*(k+1)*z.*gamma+z.*(z.*(k-1)-2*gamma*(k+1)).*coth(sqrt(1+k)*z))./(sqrt(1+k)*z.*(z+2*gamma));
% 
% dRdt = (9*k-1)/(9*sqrt(1+k))+1./z-10*coth(z*sqrt(1+k))/(9*sqrt(1+k));
% 
% figure(2)
% plot(z,dRdt)
% hold on
% plot(y,zeros(1,length(y)))


if qsave == 1

    y0 = Rinit;

    [vt,y] = ode45(@(t,y)1/y-(k*(phi_m-1)+phi_m+(k*(phi_p-1)+phi_p)*coth(sqrt(1+k)*y))/(sqrt(1+k)*(phi_p-phi_m)),0:0.05:100,y0);
       
    vR = y;

    dlmwrite(['veusz_Rvst_k_' num2str(k) '_Rinit_' num2str(Rinit)  '.txt'],[vt vR],'delimiter',' ');
    %dlmwrite(['veusz_dRdtvsR_k_' num2str(k) '_gamma_' num2str(gamma) '_Rinit_' num2str(Rinit)  '.txt'],[z' dRdt'],'delimiter',' ');
end

