close all
clear all %#ok<CLALL>
clc

addpath(genpath('../PostProcessingReqs'))

dn = 10;

fldr = "../../Results/";
Animations = true;
AnimationsIP = false;
defscale = 10;


	meshFile = fldr+"mesh_0.hdf5";
	
	%TimeData
	TData = h5read(fldr+"TimeData.hdf5", "/TimeData");
	t_Data = TData(:,1);
	dt_Data = TData(:,2);
	L_Data = TData(:,4);
	fgr = figure(99);
	subplot(2,1,1)
	t2 = repelem(t_Data,2); t2(1)=[];t2(end)=[];
	semilogy(t2, repelem(dt_Data(2:end),2));
	xlabel("$t\;[\mathrm{s}]$",'Interpreter','latex');
	ylabel("$\Delta t\;[\mathrm{s}]$", 'Interpreter','latex');
	ax = gca;
	lims = [0.0, max(t_Data)];
	ax.XLim = lims;
	subplot(2,1,2)
	plot(t_Data, L_Data)
	xlabel("$t\;[\mathrm{s}]$",'Interpreter','latex');
	ylabel("$\sum L_{\mathrm{cracks}} \;[\mathrm{m}]$", 'Interpreter','latex');
	ax = gca;
	ax.XLim = lims;
	saveFigNow(fgr, "TimeLFrac", 8, false, false, 0)


	nmax= length(t_Data);
	ind = 0;
	for i=dn:dn:nmax
		ind = ind+1;
		dataFile{ind} = fldr+"results_"+string(i)+".hdf5"; 
	end

	if (Animations)
		DataNamesNodes = {"ux","uy","Sw","phase","pw"};
		GroupName = "/"+{"3_internal","3_internal","2_internal","2_internal","2_internal","2_internal"}+"/";

		% load mesh
		XMesh = h5read(meshFile,GroupName{1}+"X"); %#ok<*UNRCH>
		YMesh = h5read(meshFile,GroupName{1}+"Y");
		
		for i=1:length(DataNamesNodes)
			fg{i} = figure(i);
			fg{i}.Position = [-40, 5, 30, 20];
			vidfile{i} = VideoWriter("Animations/"+DataNamesNodes{i}+".mp4",'MPEG-4');
			vidfile{i}.FrameRate = 10;
			open(vidfile{i});
			fg{i}.WindowState="maximized";	
		end
	
		for j=1:length(dataFile)
			for i=1:length(DataNamesNodes)
				t = h5read(dataFile{j},"/time");
		
				Data = h5read(dataFile{j},GroupName{i}+DataNamesNodes{i});
				D{i} = Data;
			end
		
			for i=1:length(DataNamesNodes)
				figure(fg{i})
				ax = gca;
				clf
				set(gcf,'color','w');
		
				if (i==1)
					PlotNodeData(XMesh, YMesh, D{i}, D{1}, D{2}, defscale);
				else
					PlotNodeData(XMesh, YMesh, D{i}, D{1}, D{2}, 0);	
				end
				if(i==3 || i==4)
					clim([0 1])
				end
				title(["t= "+string(t)+" s";DataNamesNodes{i}]);
		
				axis image;
				%camlight;
				%lighting gouraud;
				%view(3)
		
				colormap(jet)
				drawnow();
				frame = getframe(gcf);
				writeVideo(vidfile{i}, frame);
			end
		end
		
		for i=1:length(DataNamesNodes)
			close(vidfile{i});
			close(fg{i});
		end
	end


	if AnimationsIP
		DataNamesIP = {"s1","s2","s3","dam"};
		GroupName = "/"+{"2_internal","2_internal","2_internal","2_internal"}+"/";

		% load mesh
		XIPMesh = h5read(meshFile,GroupName{1}+"Xip"); %#ok<*UNRCH>
		YIPMesh = h5read(meshFile,GroupName{1}+"Yip");

		fg{1} = figure(i);
		fg{1}.Position = [-40, 5, 30, 20];
		vidfile{1} = VideoWriter("Animations/"+"Stresses.mp4",'MPEG-4');

		vidfile{1}.FrameRate = 10;
		open(vidfile{1});
		fg{1}.WindowState="maximized";	
	
		for j=1:length(dataFile)
			for i=1:length(DataNamesIP)
				t = h5read(dataFile{j},"/time");
		
				Data = h5read(dataFile{j},GroupName{i}+DataNamesIP{i});
				D{i} = Data;
			end
		
			figure(fg{1})
			clf;
			tiledlayout(3,2);
			set(gcf,'color','w');
			for i=1:length(DataNamesIP)+2
				nexttile
				
				if (i<=length(DataNamesIP))
					PlotIPData(XIPMesh, YIPMesh, D{i});
					title(["t= "+string(t)+" s";DataNamesIP{i}]);
				elseif (i==length(DataNamesIP)+1)
					PlotIPData(XIPMesh, YIPMesh, 1/sqrt(2)*sqrt((D{1}-D{2}).^2+(D{2}-D{3}).^2+(D{3}-D{1}).^2)  );
					title(["t= "+string(t)+" s";"\tau"]);
				else
					PlotIPData(XIPMesh, YIPMesh,  (D{1}+D{2}+D{3})/3  );
					title(["t= "+string(t)+" s";"p"]);
				end
		
				axis image;
				%camlight;
				%lighting gouraud;
				view(2)
		
				colormap(jet)
				%clim([-2e6, 2e6]);
				cb = colorbar;
			end
			%clim([-100 100])
			%cb = colorbar;
			%cb.Layout.Tile = "east";
			

			drawnow();
			frame = getframe(gcf);
			writeVideo(vidfile{1}, frame);
		end
		close(vidfile{1});
		close(fg{1});

	end

function PlotIPData(Xip, Yip, D)  %#ok<DEFNU>
	persistent x3d y3d F xv sel nData;

	np = 400;

	dv = D(:);
	renewData = false;
	if (isempty(x3d))
		renewData = true;
	elseif (nData~=length(Xip))
		renewData = true;
	end

	if (renewData)
		nData = length(Xip);

		xv = Xip(:);
		yv = Yip(:);

		[~, sel] = unique([xv,yv],'rows');
		xv = Xip(:); xv = xv(sel);
		yv = Yip(:); yv = yv(sel);

		[x3d, y3d] = meshgrid(linspace(min(xv),max(xv),np), linspace(min(yv),max(yv),np) );
		F = scatteredInterpolant(xv,yv,dv(sel),'linear','none');
	else 
		F.Values = dv(sel);
	end

	V3d = F(x3d(:), y3d(:));
	V3d = reshape(V3d, size(x3d));
	surf(x3d, y3d, V3d,'EdgeColor','interp','FaceColor','interp');
end



function PlotNodeData(X, Y, Data, dx, dy, defScale) %#ok<DEFNU>
	%persistent order X_el Y_el Z_el
	X_el = []; order={};

	if (isempty(X_el) || defScale~=0)
		order = [1 2 4 3];
	
		X_el = zeros(size(X,2), length(order));
		Y_el = zeros(size(X,2), length(order));
		Data_el = zeros(size(X,2), length(order));
		c=0;
		for i=1:size(X,2)
			c=c+1;
			X_el(c,:) = X(order,i)+defScale*dx(order,i);
			Y_el(c,:) = Y(order,i)+defScale*dy(order,i);
			Data_el(c,:) = Data(order,i);
		end
	end

	Data_el = zeros(size(X,2), length(order));
	c=0;
	for i=1:size(X,2)
		c=c+1;
		Data_el(c,:) = Data(order,i);
	end

	patch(X_el',Y_el',Data_el','EdgeColor','interp','FaceColor','interp');
	colorbar
end

function PlotContour(X,Y,Z,data, threshold) %#ok<DEFNU>
	persistent x3d y3d z3d F srfOutside xv yv zv sel nData;

	np = 100;

	dv = data(:);
	renewData = false;
	if (isempty(x3d))
		renewData = true;
	elseif (nData~=length(X))
		renewData = true;
	end

	if (renewData)
		nData = length(X);

		xv = X(:);
		yv = Y(:);
		zv = Z(:);

		[~, sel] = unique([xv,yv,zv],'rows');
		xv = X(:); xv = xv(sel);
		yv = Y(:); yv = yv(sel);
		zv = Z(:); zv = zv(sel);

		[x3d, y3d, z3d] = meshgrid(linspace(min(xv),max(xv),np), linspace(min(yv),max(yv),np),   linspace(min(zv),max(zv),np) );
		F = scatteredInterpolant(xv,yv,zv,dv(sel),'linear','none');

		srfOutside = convhulln([xv, yv, zv]);	
	else 
		F.Values = dv(sel);
	end

	V3d = F(x3d(:), y3d(:), z3d(:));
	%V3d(isnan(V3d)) = -10.0;
	V3d = reshape(V3d, size(x3d));
	srf = isosurface(x3d,y3d,z3d,V3d,threshold);

	p = patch(srf,'FaceColor','r');
	p.EdgeColor = "none";
	p.FaceAlpha = 0.8;
	hold on

	p2 = trisurf(srfOutside,xv, yv, zv,'FaceColor','k');
	p2.EdgeColor = "none";
	p2.FaceAlpha = 0.1;
end




function saveFigNow(fg, sname, HFig, WFig, hasColorbar, cb)
	figure(fg);
	fprintf(sname+"  ")
	ax = gca;
	ax.FontSize = 8;
	fg.Units = 'centimeters';
	if (WFig==true || WFig==false)
		if (WFig)
			fg.Position = [2 2 16 HFig];
		else
			fg.Position = [2 2 8 HFig];
		end
	else
		fg.Position = [2 2 WFig HFig];
	end
	if (hasColorbar)
		cb.Position(1) = 0.85;
		cb.Position(2) = 0.1;
		cb.Position(4) = 0.65;
		if(HFig==5)
			ax.Position(1) = 0.10;
			ax.Position(2) = 0.20;
			ax.Position(3) = 0.70;
			ax.Position(4) = 0.70;
		else
			ax.Position(1) = 0.05;
			ax.Position(2) = 0.05;
			ax.Position(3) = 0.75;
			ax.Position(4) = 0.85;
		end
	end
	set(fg,'color','w');

	drawnow();
	print(fg, "Figures/"+sname+".png",'-dpng','-r1200'); fprintf(".png  ")
	print(fg, "Figures/"+sname+".jpg",'-djpeg','-r1200'); fprintf(".jpg  ")
	print(fg, "Figures/"+sname+".eps",'-depsc','-r1200'); fprintf(".eps  ")
	print(fg, "Figures/"+sname+".svg",'-dsvg','-r1200'); fprintf(".svg  ")
	print(fg, "Figures/"+sname+".emf",'-dmeta','-r1200'); fprintf(".emf\n")
end


