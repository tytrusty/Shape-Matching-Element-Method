% example.m plots an IGES CAD-object

% Compile the c-files
makeIGESmex;

% Load parameter data from IGES-file.
[ParameterData,EntityType,numEntityType,unknownEntityType,numunknownEntityType]=iges2matlab('IGESfiles/example.igs');

% Plot the IGES object
plotIGES(ParameterData, 0, 1, 10);
%FacePlot=plotIGES(ParameterData,1,1,1000,1,10,1,'r')