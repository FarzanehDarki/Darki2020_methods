clear; close all
%% parameters
h=0.5;        %adaptation strength; WTA:(h=0.5), RIV(h=2); SIM(h=6)
g=1.5;        %inhibition strength
tau=20;       %ms
tau_H=900;    %ms
tau_I=11;     %ms
f2=18;        %Hz, flicker frequency
t_total=6000; %ms
steptime=0.1; %ms
v=10;         %max of inputs to populations representing vertical grating 
              %in the right eye and horizontal grating in the left eye
%% initial conditions
n=t_total/steptime;
E_V_Right=zeros(1,n); E_H_Left=zeros(1,n); H_V_Right=zeros(1,n);
H_H_Left=zeros(1,n); I_V_Right=zeros(1,n); I_H_Left=zeros(1,n);
E_V_Right(1)=0.000001; 
%% Solve ODEs of the Wilson model
V_Right=zeros(1,n);
H_Left=zeros(1,n);
for i=1:n-1
   xx=cos(2*pi*f2*0.001*i*steptime-15000*steptime);
   V_Right(i+1)= v/(1+exp(-10*xx)); %inputs to populations representing vertical grating in the right eye
   H_Left(i+1)= v/(1+exp(-10*xx));  %inputs to populations representing horizontal grating in the left eye
   E_V_Right(i+1)=E_V_Right(i)+steptime*(1/tau)*(-E_V_Right(i)+(100*((V_Right(i)-g*I_H_Left(i))>0)*(V_Right(i)-g*I_H_Left(i))^2)/((10+H_V_Right(i))^2+((V_Right(i)-g*I_H_Left(i))>0)*(V_Right(i)-g*I_H_Left(i))^2));
   H_V_Right(i+1)=H_V_Right(i)+steptime*(1/tau_H)*(-H_V_Right(i)+h*E_V_Right(i));
   I_V_Right(i+1)=I_V_Right(i)+steptime*(1/tau_I)*(-I_V_Right(i)+E_V_Right(i));
   E_H_Left(i+1)=E_H_Left(i)+steptime*(1/tau)*(-E_H_Left(i)+(100*((H_Left(i)-g*I_V_Right(i))>0)*(H_Left(i)-g*I_V_Right(i))^2)/((10+H_H_Left(i))^2+((H_Left(i)-g*I_V_Right(i))>0)*(H_Left(i)-g*I_V_Right(i))^2));
   H_H_Left(i+1)=H_H_Left(i)+steptime*(1/tau_H)*(-H_H_Left(i)+h*E_H_Left(i));
   I_H_Left(i+1)=I_H_Left(i)+steptime*(1/tau_I)*(-I_H_Left(i)+E_H_Left(i));
end
%% Plot time histories
t=(steptime:steptime:t_total)/1000; % convert to [s] from [ms]
figure
plot(t,E_V_Right,'b')
hold on
plot(t,E_H_Left,'b--')
title(['Flicker only case,  [J_{HL}]_{max}=[J_{VR}]_{max}=',num2str(v),',  g=',num2str(g),',  h=',num2str(h)])
ylabel('E_1 or E_2')
xlabel('Time [Sec]')