%% t-splines unequal order mesh generation
warning('off','MATLAB:imagesci:hdf5dataset:datatypeOutOfRange')

close all
addpath(genpath('../../../../../Libraries/igafem'))
clearvars
clc
	delete(gcp('nocreate'))
	parpool('threads')


cx = 14;
cy = 4;
for rat=[5]
	for dx=[1]
		for H = [200]
		
			LX = rat*H;
			LY = H;
			
			plotlim = 9999;
			layers = 0;
			Nx = ceil(LX/dx);
			Ny = ceil(LY/dx);
			
			savename = "IceCliff_L"+string(LX)+"m_H"+string(H)+"m_dx"+string(dx)+"m_"+string(cx*cy)+"c";
			disp(savename);
			surfaceMap = [0, LX; 
	  					LY, LY];
			
			mesh1 = create_anchors_simple(2, 1, Nx, Ny, surfaceMap);
			mesh1 = process_anchorsv3(mesh1);
			mesh_total = mesh1;
			
			[Nodes, NodeGroup, ElementGroup, coreData] = partition(mesh_total, cx, cy);
			saveToHDF(savename, Nodes, NodeGroup, ElementGroup, coreData)
			
			PrintInfo(mesh_total);
		end
	end
end

function PrintInfo(mesh)
	fprintf("Mesh Information\n");

	fprintf("\t Total number of Nodes: "+string(size(mesh.anchors,2))+"\n");
	for i=1:length(mesh.nodegroups)
		fprintf("\t\t "+mesh.nodegroups{i}.name+": ");
		fprintf(string(length(mesh.nodegroups{i}.nodes))+"\n");
	end

	fprintf("\t Number of Elements: "+"\n");
	for i=1:length(mesh.elementgroups)
		fprintf("\t\t "+mesh.elementgroups{i}.name+": ");
		fprintf(string(size(mesh.elementgroups{i}.elements,2))+"  ("+mesh.elementgroups{i}.type+")\n");
	end

end

%% related functions
function plot_phys_domain(mesh, offset)

    % plot cp's
    if (true)
        for an = 1:length(mesh.anchors)  
            cp(an,:) = mesh.anchors{an}.cp;        
            %text(cp(1),cp(2),int2str(an+offset),'FontSize',8)      
        end
        plot(cp(:,1),cp(:,2),'k*')
        hold on
    end
    
    if (true)
        % plot elements
        rp = 0;
        for el=1:length(mesh.elementgroups{1}.elements)
            crn = mesh.elementgroups{1}.elements{el}.corners_phys;

            x = [crn(1,1) crn(1,2) crn(1,2) crn(1,1)];
            y = [crn(2,1) crn(2,1) crn(2,2) crn(2,2)];

            rp = rp+1;
            x_rec(rp,1:4) = x;
            y_rec(rp,1:4) = y;

            xt(rp) = x(1);
            yt(rp) =(y(1)+y(3))/2;
            txt{rp}= "E"+string(el);

            %text(x(1),(y(1)+y(3))/2,"E"+string(el),'Color','blue','FontSize',10)   
        end
        patch('XData',x_rec', 'YData',y_rec','FaceColor','None','FaceAlpha',0.2,'EdgeColor','k')
        %text(xt, yt, txt, 'Color','blue','FontSize',10) 
    end
    
    title("p="+string(mesh.p)+", physical domain")
    xlabel('x [m]')
    ylabel('y [m]')
end

function vals = plot_global_basisfunc(local, xplot, yplot, plotting)

    b_spline_x = spmak(local(1,:),1);
    b_spline_y = spmak(local(2,:),1);

    Rx = fnval(b_spline_x, xplot);
    Ry = fnval(b_spline_y, yplot);
    R = (Rx'*Ry);
    if (plotting)
        [rows, cols] = find(R);
        rows = [min(rows) max(rows)];
        cols = [min(cols) max(cols)];
        surf(xplot(rows(1):rows(2)), yplot(cols(1):cols(2)), R(rows(1):rows(2), cols(1):cols(2))','FaceAlpha',0.8)
    end
    
    vals = R;
end

function saveToHDF(fname, Nodes, NodeGroup, ElementGroup, coreData)
	if exist(fname+".h5", 'file')==2
		delete(fname+".h5");
	end

	% nodes
	h5create(fname+".h5",'/nodes',size(Nodes),'Datatype','double');
	h5write(fname+".h5",'/nodes', Nodes)

	% nodegroups
	h5create(fname+".h5",'/nodegroupnames',length(NodeGroup),'Datatype','string');
	for i=1:length(NodeGroup)
		names{i} = NodeGroup{i}.name;
	end

	h5write(fname+".h5",'/nodegroupnames', string(names))
	for i=1:length(NodeGroup)
		h5create(fname+".h5",'/nodegroups/'+NodeGroup{i}.name,length(NodeGroup{i}.nodes),'Datatype','uint64');
		h5write(fname+".h5",'/nodegroups/'+NodeGroup{i}.name,uint64(NodeGroup{i}.nodes)-1);
	end

	clear names 
	% Elementgroups
	h5create(fname+".h5",'/elementgroupnames',length(ElementGroup),'Datatype','string');
	h5create(fname+".h5",'/elementgrouptypes',length(ElementGroup),'Datatype','string');
	for i=1:length(ElementGroup)
		names{i} = ElementGroup{i}.name;
		types{i} = ElementGroup{i}.type;
	end

	h5write(fname+".h5",'/elementgroupnames', string(names))
	h5write(fname+".h5",'/elementgrouptypes', string(types))
	for i=1:length(ElementGroup)
		ELS = zeros(length(ElementGroup{i}.elements), length(ElementGroup{i}.elements{1}.cps));
		DATA = zeros(length(ElementGroup{i}.elements), size(ElementGroup{i}.elements{1}.Bezier,1), size(ElementGroup{i}.elements{1}.Bezier,2));
		for el=1:length(ElementGroup{i}.elements)
			ELS(el,:) = uint64(ElementGroup{i}.elements{el}.cps);
			DATA(el,:,:) = ElementGroup{i}.elements{el}.Bezier;
		end
		h5create(fname+".h5",'/elementgroups/'+ElementGroup{i}.name,size(ELS),'Datatype','uint64');
		h5write(fname+".h5",'/elementgroups/'+ElementGroup{i}.name,ELS-1);

		h5create(fname+".h5",'/elementgroups/'+ElementGroup{i}.name+"_Data",size(DATA),'Datatype','double');
		h5write(fname+".h5",'/elementgroups/'+ElementGroup{i}.name+"_Data",DATA);
	end

	%h5create(fname+".h5",'/matchinggroups',size(matchingGroups),'Datatype','uint64');
	%h5write(fname+".h5",'/matchinggroups', matchingGroups-1);

	% partition data
	h5create(fname+".h5",'/NCores', 1,'Datatype','uint64');
	h5write(fname+".h5",'/NCores', uint64(coreData.NCores));

	for c=1:coreData.NCores
		h5create(fname+".h5",'/Partition/'+string(c-1)+'/noderange', 2,'Datatype','uint64');
		h5write(fname+".h5",'/Partition/'+string(c-1)+'/noderange',coreData.ToSave{c}.Noderange-1);

		if (isempty(coreData.ToSave{c}.Ghosts))
			h5create(fname+".h5",'/Partition/'+string(c-1)+'/Hasghosts', 1,'Datatype','uint8');
			h5write(fname+".h5",'/Partition/'+string(c-1)+'/Hasghosts', uint8(false));
		else
			h5create(fname+".h5",'/Partition/'+string(c-1)+'/Hasghosts', 1,'Datatype','uint8');
			h5write(fname+".h5",'/Partition/'+string(c-1)+'/Hasghosts', uint8(true));

			h5create(fname+".h5",'/Partition/'+string(c-1)+'/ghosts', length(coreData.ToSave{c}.Ghosts),'Datatype','uint64');
			h5write(fname+".h5",'/Partition/'+string(c-1)+'/ghosts',coreData.ToSave{c}.Ghosts-1);
		end


		h5create(fname+".h5",'/Partition/'+string(c-1)+'/Haselems', length(coreData.ToSave{c}.haselems),'Datatype','uint8');
		h5write(fname+".h5",'/Partition/'+string(c-1)+'/Haselems',uint8(coreData.ToSave{c}.haselems));

		h5create(fname+".h5",'/Partition/'+string(c-1)+'/elemrange', size(coreData.ToSave{c}.Elemrange),'Datatype','uint64');
		h5write(fname+".h5",'/Partition/'+string(c-1)+'/elemrange',coreData.ToSave{c}.Elemrange-1);

		h5create(fname+".h5",'/Partition/'+string(c-1)+'/hasnodegroup', length(coreData.ToSave{c}.hasnodegroup),'Datatype','uint8');
		h5write(fname+".h5",'/Partition/'+string(c-1)+'/hasnodegroup',uint8(coreData.ToSave{c}.hasnodegroup));

		h5create(fname+".h5",'/Partition/'+string(c-1)+'/nodegrouprange', size(coreData.ToSave{c}.NodeGrouprange),'Datatype','uint64');
		h5write(fname+".h5",'/Partition/'+string(c-1)+'/nodegrouprange',coreData.ToSave{c}.NodeGrouprange-1);
	end

	%h5disp(fname+".h5")
end
