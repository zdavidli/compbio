function [dydt] = Enz_progFunc(t,y )
%Function to calculate enzymatic reaction progress curves
%
%        kp1     kp2    kp3
% S + E <--> ES <--> EP --> E + P
%        km1     km2     
%
global kp1 km1 kp2 km2 kp3;

%Retrieve name of array variables
% Not necessary but makes reading equations easier
E=y(1);
S=y(2);
ES=y(3);
EP=y(4);
P=y(5);

% Differential equations
dE= -kp1*E*S + km1*ES +kp3*EP;
dS= -kp1*E*S + km1*ES;
dES= kp1*E*S - km1*ES - kp2*ES + km2*EP;
dEP= kp2*ES - km2*EP - kp3*EP;
dP = kp3*EP;

% Assign values to dydt vector
dydt(1)= dE;
dydt(2)= dS;
dydt(3)= dES;
dydt(4)= dEP;
dydt(5)= dP;
dydt = dydt';       %MATLAB wants a column vector

end

