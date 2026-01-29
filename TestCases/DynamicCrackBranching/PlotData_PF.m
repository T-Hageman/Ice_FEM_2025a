close all
clear all
clc

addpath(genpath('../PostProcessingReqs'))

meshFile = "../../Results/mesh_0.hdf5";
dataFile = "../../Results/results_4.hdf5";
DataNamesNodes = {"ux","uy","phase"};
GroupName = "/"+{"", "", "","","","","", "", "","","","","", "", "","","",""}+"internal/";%+"Interior/";
DataNamesIPs = {"sxx","syy","szz","sxy"};

defScale = 1*1e2;

individualfigures = false;

% Load mesh
%h5disp(meshFile)

% Load data
%h5disp(dataFile)
if (individualfigures == false)
	f = figure;
	t = tiledlayout('flow');
	t.Padding = "tight";
	t.TileSpacing = "tight";
	f.Units = "centimeters";
	f.Position = [1, 0, 50, 25];
end
read = false;
while read == false;
	try
		t = h5read(dataFile,"/time")
		
		for i=1:length(DataNamesNodes)
			X = h5read(meshFile,GroupName{i}+"X");
			Y = h5read(meshFile,GroupName{i}+"Y");
		
			Data = h5read(dataFile,GroupName{i}+DataNamesNodes{i});
			D{i} = Data;
		end
		for i=1:length(DataNamesIPs)
			Xip = h5read(meshFile,GroupName{i}+"Xip");
			Yip = h5read(meshFile,GroupName{i}+"Yip");
		
			DIP{i} = h5read(dataFile,GroupName{i}+DataNamesIPs{i});
		end
		read = true;
	catch ME
		pause(2)
	end
end

dx = D{1};
dy = D{2};
for i=1:length(DataNamesNodes)
	if (individualfigures)
		figure
	else
		nexttile
	end
	PlotNodeData(X, Y, D{i}, dx, dy, defScale);
	title(DataNamesNodes{i});
	%axis image
end

for i=1:length(DataNamesIPs)
	nexttile
	PlotIPData(Xip, Yip, DIP{i});
	title(DataNamesIPs{i});
	%axis image
	% if (i>1 && i<=4)
	% 	clim([-2e6 0.5e6]);
	% end
	% if (i==5)
	% 	clim([-1e6 1e6]);
	% end
end


function PlotNodeData(X, Y, Data, dx, dy, defScale)
if (size(X,1)==9)
	order = [1 2 3 6 9 8 7 4];
elseif (size(X,1)==6)
	order = [1 2 3];
else
	order = [1 2 4 3];
end
	X_el = zeros(size(X,2), length(order));
	Y_el = zeros(size(X,2), length(order));
	Data_el = zeros(size(X,2), length(order));

	for i=1:size(X,2)
		X_el(i,:) = X(order,i)+defScale*dx(order,i);
		Y_el(i,:) = Y(order,i)+defScale*dy(order,i);
		Data_el(i,:) = Data(order,i);
	end
	p = patch(X_el',Y_el',Data_el',Data_el','EdgeColor','interp','FaceColor','interp');
	colorbar
end

function PlotIPData(X, Y, Data)
	[xplot, yplot] = meshgrid(linspace(min(min(X)), max(max(X)),200), linspace(min(min(Y)), max(max(Y)),200));

	F = scatteredInterpolant(X(:), Y(:), Data(:));
	zplot = F(xplot, yplot);

	surf(xplot, yplot, zplot,'EdgeColor','none','FaceColor','interp');
	view(2)
	colorbar
end