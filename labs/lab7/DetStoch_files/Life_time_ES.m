
%%Calculate ES lifetime of tagged particle state from COPASI 
%  enzyme rxn simulation

%State of particle is in file, ES.dat
%
% This file has column 4 extracted from TaggedParticle.txt
%
% The values are binary instance of ES state or not 
%
clear all; % clears the memory of previous variables

% Specify data file
file = load('ES.dat');
data = file(:);

% Initialize variables
lastES = 0;
lastES = 0;
numES = 0;
numES = 0;
thisESlife = 0;
thisESlife = 0;
timeES = 0;
timeES = 0;

% Count consecutive instances of "one" in column data and sum their length
for i=1:length(data)
    if (data(i) == 1)
        if (lastES == 0);
            numES = numES + 1;
            lastES = 1;
        end
        
        if (numES > 0 && lastES == 1);
            thisESlife = thisESlife + 1;
            timeES(numES) = thisESlife;
        end    
    end
    if (data(i) == 0)
        thisESlife = 0;
        lastES = 0;
    end
end

% Plot Incidence versus lifetime
figure(5);
clf(5);

% Calculate histogram for integer data
cntES = histcounts(timeES,'BinMethod','integers',...
    'BinLimits',[0.1,99.5]);

% Calculate approximate lifetime
% Initial value = cntES(1)
% e = exp(1)
one_over_e = cntES(1)/exp(1);
for i=1:length(cntES)
    if one_over_e > cntES(i) 
        lifetime = i; 
        break
    end
end
% Since time intervals were every 0.1 second, divide by 10
scaled_lifetime = lifetime * 0.1;
fprintf( 'Lifetime = approx. %8.2f sec\n', scaled_lifetime )


% Plot it as bar graph
hES = bar(cntES);
hold on;
% Color bars blue
set(hES, 'FaceColor', 'b');
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
title(['ES Lifetime']);
xlabel('Time (s)');
ylabel('Count');
hold off;









