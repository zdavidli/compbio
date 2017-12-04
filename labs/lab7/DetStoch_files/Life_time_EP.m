
%%Calculate EP lifetime of tagged particle state from COPASI 
%  enzyme rxn simulation

%State of particle is in file, EP.dat
%
% This file has column 4 extracted from TaggedParticle.txt
%
% The values are binary instance of EP state or not 
%
clear all; % clears the memory of previous variables

% Specify data file
file = load('EP.dat');
data = file(:);

% Initialize variables
lastEP = 0;
lastEP = 0;
numEP = 0;
numEP = 0;
thisEPlife = 0;
thisEPlife = 0;
timeEP = 0;
timeEP = 0;

% Count consecutive instances of "one" in column data and sum their length
for i=1:length(data)
    if (data(i) == 1)
        if (lastEP == 0);
            numEP = numEP + 1;
            lastEP = 1;
        end
        
        if (numEP > 0 && lastEP == 1);
            thisEPlife = thisEPlife + 1;
            timeEP(numEP) = thisEPlife;
        end    
    end
    if (data(i) == 0)
        thisEPlife = 0;
        lastEP = 0;
    end
end

% Plot Incidence versus lifetime
figure(6);
clf(6);

% Calculate histogram for integer data
cntEP = histcounts(timeEP,'BinMethod','integers',...
    'BinLimits',[0.1,99.5]);

% Calculate approximate lifetime
% Initial value = cntEP(1)
% e = exp(1)
one_over_e = cntEP(1)/exp(1);
for i=1:length(cntEP)
    if one_over_e > cntEP(i) 
        lifetime = i; 
        break
    end
end
% Since time intervals were every 0.1 second, divide by 10
scaled_lifetime = lifetime * 0.1;
fprintf( 'Lifetime = approx. %8.2f sec\n', scaled_lifetime )

% Plot it as bar graph
hEP = bar(cntEP);
hold on;
% Color bars cyan
set(hEP, 'FaceColor', 'c');
xlim([0.1 100.5]);
% Put tick mark every 10th value
set(gca,'XTick',[10 20 30 40 50 60 70 80 90 100]);
% Since values were every 0.1 second, label appropriately
set(gca,'XTickLabel',[1 2 3 4 5 6 7 8 9 10]);
% Plot dashed lines to indicate 1/e value and lifetime
plot([0 lifetime], [one_over_e one_over_e],'r--',...
    [lifetime lifetime],[0 one_over_e],'r--', 'LineWidth',2);
% Axis labels, etc.
set(gca,'FontSize',15,'FontWeight','bold');
title(['EP Lifetime']);
xlabel('Time (s)');
ylabel('Count');
hold off;









