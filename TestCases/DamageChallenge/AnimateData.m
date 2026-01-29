close all
clear all
clc

addpath(genpath('../IceFEM/TestCases/PostProcessingReqs'))

DoVids = true;
DoLines = true;


nmax= 100;
dn = 1;

dscale = 50;

i=0;
for split=["SpecStrains","VolStrains"]
    for mesh=["quadratic_025mm", "cubic_025mm", "triangleLinear_025mm", "triangleQuadratic_025mm", "quadratic_0125mm", "cubic_0125mm", "triangleLinear_0125mm", "triangleQuadratic_0125mm"]
		for dist=["Linear","Quadratic","HO_Linear","HO_Quadratic"]
            for cse=["I","II"]
                if ((dist=="HO_Linear" || dist=="HO_Quadratic") && (mesh=="triangleLinear_025mm" || mesh=="triangleQuadratic_025mm" || mesh=="triangleLinear_0125mm" || mesh=="triangleQuadratic_0125mm"))
                    
				else
					i = i+1;
					casenames{i} = cse+"_"+split+"_"+dist+"_"+mesh;
				end
			end
		end
	end
end

caseMax = 2; %length(casenames)

clrs = distinguishable_colors(caseMax);

if DoLines
	figure(451)
	subplot(2,1,1)
	for z=1:caseMax
		casename = casenames{z};
		casename_filtered = replace(casename,"_"," ");
		tfile = "./Results/"+casename+"/Results/TimeData.hdf5";

		t = h5read(tfile,"/TimeData");
		u_y = t(:,1)*1e3;
		F_y = -t(:,6)*12.7e-3;
		
		if rem(z, 2) == 1
			subplot(2,1,1)
			plot(u_y, F_y,'Color',clrs((z+1)/2,:), 'DisplayName',casename_filtered);
		else
			subplot(2,1,2)
			plot(u_y, F_y,'Color',clrs(z/2,:), 'DisplayName',casename_filtered);
		end
		hold on
	end

	xlabel('u_y [mm]')
	ylabel('F [N]')
	legend('Location','eastoutside')
end

if DoVids
	for z=1:caseMax
		casename = casenames{z}
		fldr = "./Results/"+casename+"/Results/";
		
		meshFile = fldr+"mesh_0.hdf5";

		ind = 0;
		for i=1:1:nmax
			ind = ind+1;
			dataFile{ind} = fldr+"results_"+string(i)+".hdf5";
		end
		DataNamesNodes = {"ux","uy","phase"};
		GroupName = "/"+{"internal","internal","internal"}+"/";
		
		fg = figure(1);
			fg.Units = "centimeters";
			fg.Position = [-40, 5, 30, 20];
			vidfile = VideoWriter(casename+".mp4",'MPEG-4');
			vidfile.FrameRate = 10;
			open(vidfile);
			fg.WindowState="maximized";
		
		% Load data
		for j=dn:dn:nmax
			figure(fg)
			clf
			tl = tiledlayout(2,1);
		
			t = h5read(dataFile{j},"/time");
			disp(t)
			title(tl,["t= "+string(t)+" s";""]);
		
			for i=1:length(DataNamesNodes)
				Data = h5read(dataFile{j},GroupName{i}+DataNamesNodes{i});
				D{i} = Data;
				X = h5read(meshFile,GroupName{i}+"X");
				Y = h5read(meshFile,GroupName{i}+"Y");
			
				if (i>1.5)
					if (i==2)
						X = X+D{1}*dscale;
						Y = Y+D{2}*dscale;
						D{2} = D{2}*1000;
					end
					nexttile
					PlotNodeData(X, Y, D{i});
				end
				title(DataNamesNodes{i});
				axis image
				colormap(jet)
				if (i==2)
					xlim([-0.01 0.09])
					ylim([-20e-3 40e-3])
				end
				if (i==3)
					clim([0 1])
				end
			end
			
			drawnow();
			pause(1);
			drawnow();
			frame = getframe(gcf);
			writeVideo(vidfile, frame);
		end
		close(vidfile);
		close(fg);
	end
end


function PlotNodeData(X, Y, Data)
	if (size(X,1)==4) %surface
		order = [1 2 4 3];
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
