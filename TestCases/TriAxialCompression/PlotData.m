close all
clear all
clc

addpath(genpath('../PostProcessingReqs'))

meshFile = "../../Results/mesh_0.hdf5";
dataFile = "../../Results/results_60.hdf5";
DataNamesNodes = {"ux","uy","uz","phase","T"};
GroupName = "/"+{"", "", "","","","","", "", "","","","","", "", "","","",""}+"internal/";

defScale = 1e2;

individualfigures = false;

% Load data
if (individualfigures == false)
	f = figure;
	t = tiledlayout('flow');
	t.Padding = "tight";
	t.TileSpacing = "tight";
	f.Units = "centimeters";
	f.Position = [1, 0, 50, 25];
end

t = h5read(dataFile,"/time")

for i=1:length(DataNamesNodes)
	X = h5read(meshFile,GroupName{i}+"X");
	Y = h5read(meshFile,GroupName{i}+"Y");
	Z = h5read(meshFile,GroupName{i}+"Z");

	Data = h5read(dataFile,GroupName{i}+DataNamesNodes{i});
	D{i} = Data;
end
read = true;

dx = D{1};
dy = D{2};
dz = D{3};
for i=1:length(DataNamesNodes)+1
	if (individualfigures)
		figure
	else
		nexttile
	end
	if (i<=5)
		PlotNodeData(X, Y, Z, D{i}, dx, dy, dz, defScale);
		title(DataNamesNodes{i});
	elseif (i==6)
		PlotContour(X, Y, Z, D{4}, 0.5);
	end
	%axis image
end




DataNamesIP = {"exxp","eyyp","ezzp","exyp","eyzp","exzp"};
GroupName = "/"+{"internal","internal","internal","internal","internal","internal","internal"}+"/";

% load mesh
XIPMesh = h5read(meshFile,GroupName{1}+"Xip"); %#ok<*UNRCH>
YIPMesh = h5read(meshFile,GroupName{1}+"Yip");
ZIPMesh = h5read(meshFile,GroupName{1}+"Zip");

for i=1:length(DataNamesIP)
	t = h5read(dataFile,"/time");
	Data = h5read(dataFile,GroupName{i}+DataNamesIP{i});
	D{i} = Data;
end

figure()
clf;
tiledlayout(2,4);
set(gcf,'color','w');
for i=1:length(DataNamesIP)
	nexttile

	PlotIPDataSlice(XIPMesh, YIPMesh, ZIPMesh, D{i},"xz");
	title(["t= "+string(t)+" s";DataNamesIP{i}]);

	axis image;
	%camlight;
	%lighting gouraud;
	view(2)

	colormap(jet)
	%clim([-2e6, 2e6]);
	cb = colorbar;
end

figure()
clf;
tiledlayout(2,4);
set(gcf,'color','w');
for i=1:length(DataNamesIP)
	nexttile

	PlotIPDataSlice(XIPMesh, YIPMesh, ZIPMesh, D{i},"yz");
	title(["t= "+string(t)+" s";DataNamesIP{i}]);

	axis image;
	%camlight;
	%lighting gouraud;
	view(2)

	colormap(jet)
	%clim([-2e6, 2e6]);
	cb = colorbar;
end











function PlotNodeData(X, Y, Z, Data, dx, dy, dz, defScale) %#ok<DEFNU>
	%persistent order X_el Y_el Z_el
	X_el = []; order={};

	if (isempty(X_el) || defScale==0)
		order{1} = [1 2 4 3];
		order{2} = [5 6 8 7];
		order{3} = [1 2 6 5];
		order{4} = [2 4 8 6];
		order{5} = [3 4 8 7];
		order{6} = [1 3 7 5];
	
		X_el = zeros(size(X,2)*length(order), length(order{1}));
		Y_el = zeros(size(X,2)*length(order), length(order{1}));
		Z_el = zeros(size(X,2)*length(order), length(order{1}));
		Data_el = zeros(size(X,2)*length(order), length(order{1}));
		c=0;
		for i=1:size(X,2)
			for j=1:length(order)
				c=c+1;
				X_el(c,:) = X(order{j},i)+defScale*dx(order{j},i);
				Y_el(c,:) = Y(order{j},i)+defScale*dy(order{j},i);
				Z_el(c,:) = Z(order{j},i)+defScale*dz(order{j},i);
				Data_el(c,:) = Data(order{j},i);
			end
		end
	end

	Data_el = zeros(size(X,2)*length(order), length(order{1}));
	c=0;
	for i=1:size(X,2)
		for j=1:length(order)
			c=c+1;
			Data_el(c,:) = Data(order{j},i);
		end
	end

	patch(X_el',Y_el',Z_el',Data_el','EdgeColor','interp','FaceColor','interp');
	colorbar
end

function PlotContour(X,Y,Z,data, threshold)
	np = 300;
	xv = X(:);
	yv = Y(:);
	zv = Z(:);
	dv = data(:);

	[x3d, y3d, z3d] = meshgrid(linspace(min(xv),max(xv),np), linspace(min(yv),max(yv),np),   linspace(min(zv),max(zv),np) );
	F = scatteredInterpolant(xv,yv,zv,dv,'linear','none');
	V3d = F(x3d(:), y3d(:), z3d(:));
	V3d = reshape(V3d, size(x3d));
	srf = isosurface(x3d,y3d,z3d,V3d,threshold);
	srfOutside = convhulln([xv, yv, zv]);

	p = patch(srf,'FaceColor','r');
	p.EdgeColor = "none";
	p.FaceAlpha = 0.8;
	hold on

	p2 = trisurf(srfOutside,xv, yv, zv,'FaceColor','k');
	p2.EdgeColor = "none";
	p2.FaceAlpha = 0.1;
end

function PlotIPDataSlice(Xip, Yip, Zip, D, slic)  %#ok<DEFNU>
	persistent x3d y3d z3d F xv zv sel nData;

	np = 100;

	dv = D(:);
	renewData = true;
	if (isempty(x3d))
		renewData = true;
	elseif (nData~=length(Xip))
		renewData = true;
	end

	if (renewData)
		nData = length(Xip);

		xv = Xip(:);
		yv = Yip(:);
		zv = Zip(:);

		[~, sel] = unique([xv,yv,zv],'rows');
		xv = Xip(:); xv = xv(sel);
		yv = Yip(:); yv = yv(sel);
		zv = Zip(:); zv = zv(sel);

		if (slic=="xz")
			[x3d, z3d] = meshgrid(linspace(min(xv),max(xv),np), linspace(min(zv),max(zv),np) );
			y3d = x3d*0;
		else
			[y3d, z3d] = meshgrid(linspace(min(yv),max(yv),np), linspace(min(zv),max(zv),np) );
			x3d = y3d*0;
		end
		F = scatteredInterpolant(xv,yv,zv,dv(sel),'linear','none');
	else 
		F.Values = dv(sel);
	end

	V3d = F(x3d(:), y3d(:), z3d(:));
	V3d = reshape(V3d, size(x3d));
	if (slic=="xz")
		surf(x3d, z3d, V3d,'EdgeColor','interp','FaceColor','interp');
	else
		surf(y3d, z3d, V3d,'EdgeColor','interp','FaceColor','interp');
	end
end