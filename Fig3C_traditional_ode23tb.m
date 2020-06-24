clear; close all
%% parameters
h=1;           %adaptation strength; WTA:(h=1), RIV(h=4.3); SIM(h=15)
g=1.5;         %inhibition strength
tau=20;        %ms
tau_H=900;     %ms
tau_I=11;      %ms
t_total=20000; %simulation time in [ms]
J_HL=10;       %inputs to populations representing horizontal grating in the left eye
J_VR=10;       %inputs to populations representing vertical grating in the right eye
%% Wilson model 
Wilson=@(t,u)[...
    (1/tau)*(-u(1,:)+(100*((J_HL-g*u(6,:))>0)*(J_HL-g*u(6,:))^2)/((10+u(2,:))^2+((J_HL-g*u(6,:))>0)*(J_HL-g*u(6,:))^2));
    (1/tau_H)*(-u(2,:)+h*u(1,:));
    (1/tau_I)*(-u(3,:)+u(1,:));
    (1/tau)*(-u(4,:)+(100*((J_VR-g*u(3,:))>0)*(J_VR-g*u(3,:))^2)/((10+u(5,:))^2+((J_VR-g*u(3,:))>0)*(J_VR-g*u(3,:))^2));
    (1/tau_H)*(-u(5,:)+h*u(4,:));
    (1/tau_I)*(-u(6,:)+u(4,:))];
%% Solve ODEs with ode23tb
x0=zeros(6,1); x0(1,1)=0.000001; %initial conditions
tspan=[0 t_total];
options = odeset('AbsTol',1e-7,'RelTol',1e-7);
[t,xt] = ode23tb(Wilson, tspan, x0, options);
t=t/1000; % convert to [s] from [ms]
%% Plot time histories
figure
plot(t,xt(:,1),'b')
hold on
plot(t,xt(:,4),'b--')
title(['Traditional case,  J_{HL}=J_{VR}=',num2str(J_HL),',  g=',num2str(g),',  h=',num2str(h)])
ylabel('E_1 or E_2')
xlabel('Time [Sec]')