close all
clear all
clc

addpath(genpath('../PostProcessingReqs'))

dn = 10;

fldrs = {"../.."};
p = [5e6];


for zzz=1:length(fldrs)
	fileNamePrefix{zzz} = replace(fldrs{zzz},"/",""); %#ok<*SAGROW>
	DispNameID{zzz} = replace(fileNamePrefix{zzz},"_"," ");
end

MaxResults = 1%length(fldrs);

Animations = true;
AnimationsIP = false;
ContourOverLay = false;
defscale = 100;

R = 10e-2/2;

m.Young = 5e9;
m.Poisson = 0.35;
m.Gc = 4000e-3;
m.l = 10.0e-3;
m.ereff = 5.0e-4;
m.creduce = 2e-6;

m.i = [1 1 1 0];
m.J2Mat = eye(4) - 1.0/3.0*m.i'*m.i;J2Mat(4,4)=2;
m.Shear= m.Young/(2.0*(1+m.Poisson));
m.Lame = m.Young*m.Poisson/((1.0+m.Poisson)*(1.0-2.0*m.Poisson));
m.Bulk = m.Young/(3*(1-2*m.Poisson));
m.D = m.Lame*(m.i'*m.i)+2*m.Shear*eye(4);
m.D(4,4) = m.Shear;
m.Dinv = inv(m.D);
m.ab = 1/4;
m.dDam = 2;
m.contourcrit = 0.5 * m.ab * m.dDam * m.Gc/m.l;
















ca = distinguishable_colors(20);

for zzz=1:MaxResults
	f = fldrs{zzz};
	fldr = f+"/Results/";
	meshFile = fldr+"mesh_0.hdf5";
	
	%TimeData
	TData = h5read(fldr+"TimeData.hdf5", "/TimeData");
	U_Ext = TData(:,7);
	P_Data = TData(:,6);
	figure(99)
	plot(U_Ext*1e3, P_Data*1e-6,'Color',ca(zzz,:),'DisplayName',DispNameID{zzz});
	xlabel("$U_\mathrm{ext}\;[\mathrm{mm}]$",'Interpreter','latex');
	ylabel("$p_\mathrm{ax}\;[\mathrm{MPa}]$", 'Interpreter','latex');
	hold on
	p_ax = [min(P_Data(end-50:end)) max(P_Data(end-50:end))];
	[~, iBreak] = min(P_Data);

	[~,i] = max(abs(p_ax-p(zzz)));
	pMax(zzz) = p_ax(i);


	ax = gca;
	xRange = ax.XLim;
	hold on
	plot(xRange, [pMax(zzz) pMax(zzz)]*1e-6, '-.','Color',ca(zzz,:),'HandleVisibility','off')
	legend('Location','southoutside');

	nmax= length(U_Ext);
	ind = 0;
	for i=dn:dn:nmax
		ind = ind+1;
		dataFile{ind} = fldr+"results_"+string(i)+".hdf5"; 
	end

	if (Animations)
		DataNamesNodes = {"ux","uy","uz","phase"};
		GroupName = "/"+{"internal","internal","internal","internal"}+"/";

		% load mesh
		XMesh = h5read(meshFile,GroupName{1}+"X"); %#ok<*UNRCH>
		YMesh = h5read(meshFile,GroupName{1}+"Y");
		ZMesh = h5read(meshFile,GroupName{1}+"Z");
		
		for i=1:5
			fg{i} = figure(i);
			fg{i}.Position = [-40, 5, 30, 20];
			if (i<=4)
				vidfile{i} = VideoWriter("Animations/"+fileNamePrefix{zzz}+"_"+DataNamesNodes{i}+".mp4",'MPEG-4');
			else
				vidfile{i} = VideoWriter("Animations/"+fileNamePrefix{zzz}+"_Contour.mp4",'MPEG-4');
			end
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
		
			for i=1:5
				figure(fg{i})
				ax = gca;
				clf
				set(gcf,'color','w');
		
				if (i<=3)
					PlotNodeData(XMesh, YMesh, ZMesh, D{i}, D{1}, D{2}, D{3}, defscale);
					title(["t= "+string(t)+" s";DataNamesNodes{i}]);
				elseif (i==4)
					PlotNodeData(XMesh, YMesh, ZMesh, D{i}, D{1}, D{2}, D{3}, 0);	
				else
					PlotContour(XMesh,YMesh,ZMesh,D{4}, 0.5)
					title("t= "+string(t)+" s");
				end
		
				axis image;
				%camlight;
				%lighting gouraud;
				view(3)
		
				colormap(jet)
				drawnow();
				frame = getframe(gcf);
				writeVideo(vidfile{i}, frame);
			end
		end
		
		for i=1:5
			close(vidfile{i});
			close(fg{i});
		end
	end


	if AnimationsIP
		DataNamesIP = {"exxp","eyyp","ezzp","exyp","eyzp","exzp"};
		GroupName = "/"+{"internal","internal","internal","internal","internal","internal","internal"}+"/";

		% load mesh
		XIPMesh = h5read(meshFile,GroupName{1}+"Xip"); %#ok<*UNRCH>
		YIPMesh = h5read(meshFile,GroupName{1}+"Yip");
		ZIPMesh = h5read(meshFile,GroupName{1}+"Zip");

		fg{1} = figure(i);
		fg{1}.Position = [-40, 5, 30, 20];
		vidfile{1} = VideoWriter("Animations/"+fileNamePrefix{zzz}+"_Stresses.mp4",'MPEG-4');

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
			tiledlayout(2,4);
			set(gcf,'color','w');
			for i=1:length(DataNamesIP)
				nexttile

				PlotIPDataSlice(XIPMesh, YIPMesh, ZIPMesh, D{i});
				title(["t= "+string(t)+" s";DataNamesIP{i}]);
		
				axis image;
				%camlight;
				%lighting gouraud;
				view(2)
		
				colormap(jet)
				%clim([-2e6, 2e6]);
				cb = colorbar;
			end
			clim([-100 100])
			%cb = colorbar;
			%cb.Layout.Tile = "east";
			

			drawnow();
			frame = getframe(gcf);
			writeVideo(vidfile{1}, frame);
		end
		close(vidfile{1});
		close(fg{1});

	end
end

fg = figure(99);
saveFigNow(fg, "ForceDisplacement", 8, false, false, 0)

if (ContourOverLay)
	fg = figure(98);
	plot(p(1:MaxResults)*1e-6, (pMax(1:MaxResults)-p(1:MaxResults))*1e-6, 'k*','DisplayName','Simulations')
	hold on

	nsteps = 201;	
	prange = [-1.0e7 1.0e7];
	[rad,ax] = meshgrid(linspace(prange(1),5e6,nsteps),linspace(prange(1),prange(2),nsteps));

	for j=1:nsteps
		for k=1:nsteps
			stressMat(:,j,k) = [rad(j,k);rad(j,k);ax(j,k)+rad(j,k); 0];
			strainMat(:,j,k) = m.Dinv*stressMat(:,j,k);
		end
	end

	E = zeros([nsteps,nsteps]);
	for j=1:nsteps
		for k=1:nsteps
				stress = stressMat(:,j,k);
				strain = strainMat(:,j,k);
				E(j,k) = IceSplit(m, strain);
		end
	end
	M = contourf(rad*1e-6,ax*1e-6,E,[m.contourcrit 100000*m.contourcrit],'k','FaceAlpha',0.25,'FaceColor','k','HandleVisibility','off');
	hold on
	plot(NaN,NaN,'k','DisplayName',"Analytic");

	xlabel('$p_\mathrm{rad}\;[\mathrm{MPa}]$','Interpreter','latex')
	ylabel('$p_\mathrm{ax}\;[\mathrm{MPa}]$','Interpreter','latex')
	legend('Location','southoutside','NumColumns',2)
	grid minor
	saveFigNow(fg, "StrengthSurface", 8, false, false, 0)
end



function PlotIPDataSlice(Xip, Yip, Zip, D)  %#ok<DEFNU>
	persistent x3d y3d z3d F xv zv sel nData;

	np = 100;

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
		zv = Zip(:);

		[~, sel] = unique([xv,yv,zv],'rows');
		xv = Xip(:); xv = xv(sel);
		yv = Yip(:); yv = yv(sel);
		zv = Zip(:); zv = zv(sel);

		[x3d, z3d] = meshgrid(linspace(min(xv),max(xv),np), linspace(min(zv),max(zv),np) );
		y3d = x3d*0;
		F = scatteredInterpolant(xv,yv,zv,dv(sel),'linear','none');
	else 
		F.Values = dv(sel);
	end

	V3d = F(x3d(:), y3d(:), z3d(:));
	V3d = reshape(V3d, size(x3d));
	surf(x3d, z3d, V3d,'EdgeColor','interp','FaceColor','interp');
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


function G = IceSplit(m, strain)
	evol = m.i*strain;
	edev = strain-1/3*evol*m.i';

	G = 0.5*m.Bulk*max(0.0,evol)^2 + m.Shear * edev'*[1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 0.25]*edev;
	G = G - 0.5*m.Bulk*m.creduce*(min(0.0,evol))^2/(m.ereff^2+(min(0.0,evol))^2);
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


