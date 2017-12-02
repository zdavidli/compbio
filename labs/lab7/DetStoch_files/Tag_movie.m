
%%Plots state of tagged particle from COPASI enzyme rxn simulation

%State of particle is in file, TaggedParticle.txt
%
% This file has five columns
%
%Time   S    ES    EP    P
% 0	    1	 0	   0	 0
% 0.5	1	 0	   0	 0
% 1	    0	 1	   0	 0
% 1.5	0	 0	   1	 0
% 2	    0	 0	   0	 1
% 2.5   0	 0	   0	 1
%etc
%
% Time is in seconds, each column is binary instance of that 
%    state of particle
%
clear all; % clears the memory of previous variables

% Specify data file
filename='TaggedParticle.txt';
delimiterIn = '\t';     % Columns separated by tabs
headerlinesIn = 1';     % First row is labels

% Read in data from file
A=importdata(filename,delimiterIn,headerlinesIn);

%Check column labels
disp(A.colheaders)

% Transform data in binary instance for each of S, ES, EP, P states
for i=1:length(A.data(:,1))
    x(i) = A.data(i,1);
    if (A.data(i,2) == 1)
        y(i) = 1;
    end
    if (A.data(i,3) == 1)
        y(i) = 2;
    end
    if (A.data(i,4) == 1)
        y(i) = 3;
    end
    if (A.data(i,5) == 1)
        y(i) = 4;
    end
end

% Plot the evolution of states with time
figure(3)
clf(3)
xlim([min(x(:)) max(x(:))])
ylim([min(y(:)) max(y(:))])

% Trick to instantiate plot - plot everything at once with invisible dots
plot(x,y,'w.')
hold on

% Label fonts, etc.
states = ['Subs';'EnzS';'EnzP';'Prod'];    %State labels for y axis
set(gca,'FontSize',15,'FontWeight','bold')
set(gca,'XTick',[0:1000:5000])    % Hard wired for current COPASI file
set(gca,'YTick',[1:4], 'YTickLabel', states)    % Four possible states
title([ 'Tagged Sustrate'],'FontSize',15','FontWeight','bold');
xlabel('Time (s)')

% Now animate in binary existance for each state as circles
for k = 1:length(A.data(:,1))
    if y(k) ~= 1
       plot(x(k),y(k),'ro');
       pause(0.05);     % Necessary but made as short as effective
    end
    
    % Stop when product (P) is formed
    if y(k) > 3
        break;
    end
end

% Add lines to plot
plot(x,y,'b-','LineWidth',2);
hold off






