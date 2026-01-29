close all
clear all
clc

addpath(genpath('../PostProcessingReqs'))

fldr = "../../Results/"

meshFile = fldr+"mesh_0.hdf5";
nmax= 440;
dn = 1;

defscale = 1000;

ind = 0;
for i=1:1:nmax
	ind = ind+1;
	dataFile{ind} = fldr+"results_"+string(i)+".hdf5";
end
DataNamesNodes = {"ux","uy","phase"};
GroupName = "/"+{"internal","internal","internal"}+"/";

individualfigures = false;
figSets = [true];

% Load mesh
%h5disp(meshFile)

if (figSets(1))
	if (individualfigures)
		for i=1:3
			fg{i} = figure(i);
			fg{i}.Position = [-40, 5, 30, 20];
			vidfile{i} = VideoWriter("IceSheet_"+DataNamesNodes{i}+".mp4",'MPEG-4');
			vidfile{i}.FrameRate = 10;
			open(vidfile{i});
			fg{i}.WindowState="maximized";		
		end
	else
		fg = figure(1);
			fg.Units = "centimeters";
			fg.Position = [-40, 5, 30, 20];
			vidfile = VideoWriter("Results.mp4",'MPEG-4');
			vidfile.FrameRate = 10;
			open(vidfile);
			fg.WindowState="maximized";
	end
end

% Load data
for j=dn:dn:nmax
	if (individualfigures==false && figSets(1))
		figure(fg)
		clf
		tl = tiledlayout('flow');
		set(gcf,'color','w');

		t = h5read(dataFile{j},"/time")
		title(tl,["t= "+string(t)+" s";""]);
	end

	for i=1:length(DataNamesNodes)
		t = h5read(dataFile{j},"/time");

		Data = h5read(dataFile{j},GroupName{i}+DataNamesNodes{i});
		D{i} = Data;
		if (i==4)
			D{i}(D{i}<-0.5e6)=-0.5e6;
		end
	end

	if(t<=0)
		XRef = D{1};
		YRef = D{2};
	end

	XInterior = h5read(meshFile,GroupName{1}+"X") + defscale*(D{1}-XRef);
	YInterior = h5read(meshFile,GroupName{1}+"Y") + defscale*(D{2}-YRef);

	XInteriorND = h5read(meshFile,GroupName{1}+"X");
	YInteriorND = h5read(meshFile,GroupName{1}+"Y");

	for i=1:3
		if (figSets(1))
			if (individualfigures == false)
				nexttile
			else
				figure(fg{i})
				clf
				set(gcf,'color','w');
			end
			if (i<2.5)
				PlotNodeData(XInterior, YInterior, D{i});
			else
				PlotNodeData(XInteriorND, YInteriorND, D{i});
				clim([0 1]);
			end

			title(DataNamesNodes{i});
			axis image
		end

		if (figSets(1))
			colormap(jet)
	
			if (individualfigures)
				title(["t= "+string(t)+" s";DataNamesNodes{i}]);
				drawnow();
				frame = getframe(gcf);
				writeVideo(vidfile{i}, frame);
			end
		end
	end
	
	if (individualfigures==false && figSets(1))
		drawnow();
		frame = getframe(gcf);
		writeVideo(vidfile, frame);
	end

end

if (figSets(1))
	if (individualfigures)
		for i=1:length(DataNamesNodes)
			close(vidfile{i});
		end
	else
		close(vidfile);
	end
end

function PlotCombined(Xint, Yint, Xsurf, Ysurf, DataArray)
	%Fluid  Phase (3) / pressure (4) / Sw (5) / Poros(6)
	order = [1 2 4 3];

	npoints = 0;
	V = zeros(size(Xint,2)*4,3);
	Vcol = zeros(size(Xint,2)*4,1);
	Vopacity = zeros(size(Xint,2)*4,1);
	F = zeros(size(Xint,2),4);
	for i=1:size(Xint,2)
		for j=1:4
			npoints = npoints+1;
			V(npoints,1) = Xint(j,i);
			V(npoints,2) = Yint(j,i);
			V(npoints,3) = 0.0;

			Vcol(npoints,1) =  DataArray{4}(j,i);
			Vopacity(npoints,1) =  max(0.0,min(1.0,DataArray{3}(j,i).^2+DataArray{6}(j,i)/0.01*DataArray{5}(j,i)));
		end

		F(i,1:4) = npoints-4+order;
	end
	p = patch('Faces',F,'Vertices',V,'FaceVertexCData',Vcol,'EdgeColor','none','FaceColor','interp');
	p.FaceVertexAlphaData = Vopacity;
	p.FaceAlpha = 'Interp';
	colorbar
	clim([0, 3.0e6])

	hold on
	
	V(:,1) = -V(:,1);
	p2 = patch('Faces',F,'Vertices',V,'FaceVertexCData',Vcol,'EdgeColor','none','FaceColor','interp');
	p2.FaceVertexAlphaData = Vopacity;
	p2.FaceAlpha = 'Interp';

	%Surface Water Height
	clear X_el Y_el
	BottomLineX = [];
	BottomLineY = [];
	for i=1:size(Xsurf,2)
		for j=2:-1:1
			BottomLineX(end+1) = Xsurf(j,i);
			BottomLineY(end+1) = Ysurf(j,i);
	
			X_el(i,j) = Xsurf(j,i);
			Y_el(i,j) = Ysurf(j,i);

			X_el(i,5-j) = Xsurf(j,i);
			Y_el(i,5-j) = Ysurf(j,i)+DataArray{7}(j,i);	
		end
	end
	patch(X_el',Y_el',0*Y_el','EdgeColor','none','FaceColor','b');
	patch(-X_el',Y_el',0*Y_el','EdgeColor','none','FaceColor','b');

	%Ice outline
	plot(BottomLineX,BottomLineY,'k');
	plot(-BottomLineX,BottomLineY,'k');
	plot(BottomLineX,0*BottomLineX,'k');
	plot(-BottomLineX,0*BottomLineX,'k');
	plot([BottomLineX(end),BottomLineX(end)],[0,BottomLineY(end)],'k');
	plot([-BottomLineX(end),-BottomLineX(end)],[0,BottomLineY(end)],'k');
end


function PlotNodeData(X, Y, Data)
	if (size(X,1)==4) %surface
		order = [1 2 4 3];
		X_el = zeros(size(X,2),4);
		Y_el = zeros(size(X,2),4);
		Data_el = zeros(size(X,2),4);
		for i=1:size(X,2)
			X_el(i,:) = X(order,i);
			Y_el(i,:) = Y(order,i);
			Data_el(i,:) = Data(order,i);
		end
		patch(X_el',Y_el',Data_el',Data_el','EdgeColor','interp','FaceColor','interp');
		colorbar
	end
	if (size(X,1)==3) %surface
		order = [1 2 3];
		for i=1:size(X,2)
			X_el(i,:) = X(order,i);
			Y_el(i,:) = Y(order,i);
			Data_el(i,:) = Data(order,i);
		end
		patch(X_el',Y_el',Data_el',Data_el','EdgeColor','interp','FaceColor','interp');
		colorbar
	end
	if (size(X,1)==2) %lines
		BottomLineX = [];
		BottomLineY = [];
		for i=1:size(X,2)
			for j=1:2
				BottomLineX(end+1) = X(j,i);
				BottomLineY(end+1) = Y(j,i);
		
				X_el(i,j) = X(j,i);
				Y_el(i,j) = Y(j,i);
				Data_el(i,j) = Data(j,i);
	
				X_el(i,5-j) = X(j,i);
				Y_el(i,5-j) = Y(j,i)+Data(j,i);
				Data_el(i,5-j) = Data(j,i);			
			end
		end
		plot(BottomLineX,BottomLineY,'k');
		hold on
		patch(X_el',Y_el',Data_el',Data_el','EdgeColor','none','FaceColor','interp');
		colorbar	
	end
end
