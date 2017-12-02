%Enzyme Reaction Progress Curve Plot
% File Enz_progress.m
% Calls Enz_progFunc.m to calculate the rate equations
% Plots S, (ES*10), (EP*10), P versus time
%
%
%        kp1     kp2    kp3
% S + E <--> ES <--> EP --> E + P
%        km1     km2     
%
%Initial parameter values (in 然 units for 2nd order)
%kp1=0.1;          然^-1 s^-1
%km1=2.0;          s^-1
%kp2=1.0;          s^-1
%km2=0.5;          s^-1
%kp3=0.04;         s^-1

clear all; % clears the memory of previous variables

%Rate and equilibrium constants  
global kp1 km1 kp2 km2 kp3;     % Won't have to pass as argument
kp1=0.1;
km1=2.0;
kp2=1.0;  
km2=0.5;
kp3=0.04;

%Initial concentrations
S0=50.0;           % 然 Substrate concentration
E0=1.0;            % 然 Enzyme concentration
ES0=0;         
EP0=0;
P0=0;           

%Time steps (seconds)
tmax=3000;
N=600;
t=linspace(0, tmax, N); %Makes vector of N increments from 0 to tmax

%Make a vector of initial concentrations to pass to function
y0=[E0;S0;ES0;EP0;P0];

%Call routine to solve ODE
% Pass time step vector and matrix of initial values
% Since at near-steady state the intermediate reaction rates are hidden
% this is considered a "stiff" differential equation problem. Hence,
% we use a different solver, ode23s, rather than ode45.
[t,y]=ode23s('Enz_progFunc',t,y0);

% Retrieve the components from the returned y matrix
E=y(:,1);   % This is a column vector = all rows in column 1
S=y(:,2);   % All rows in column 2
ES=y(:,3);  % All rows in column 3
EP=y(:,4);  % All rows in column 4
P=y(:,5);   % All rows in column 5
% Multiply complex concentrationsfor plotting  
ES=(ES)*10.0;
EP=(EP)*10.0;

%Plot results
figure(1)
clf(1)
plot(t,S,'g-',t,ES,'b-',t,EP,'c-',t,P,'r-','LineWidth',2);

%Set label size, etc.
set(gca,'FontSize',15,'FontWeight','bold')
xlabel('Time (sec)'); ylabel('Concentration (uM)');
title(['Enzymatic Reaction Progress Curves'])

%Locate the legend in best position to not overlap data
legend('Substrate','(ES)*10','(EP)*10','Product','Location','Best')
%Turn off box around legend
legend BOXOFF

