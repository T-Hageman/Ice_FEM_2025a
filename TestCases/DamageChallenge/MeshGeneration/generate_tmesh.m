%% t-splines unequal order mesh generation
close all
addpath(genpath('../../../../Libraries/igafem'))
clearvars
clc
	delete(gcp('nocreate'))
	parpool('threads')

for i=1:8
	switch i
		case 1
			type = "B";
			savename = "quadratic_025mm";
			dx = 0.25e-3;
			order = 2;
		case 2
			type = "B";
			savename = "cubic_025mm";
			dx = 0.25e-3;
			order = 3;
		case 3
			type = "Triangle";
			savename = "triangleLinear_025mm";
			dx = 0.25e-3;
			order = 1;
		case 4
			type = "Triangle";
			savename = "triangleQuadratic_025mm";
			dx = 0.25e-3;
			order = 2;
		case 5
			type = "B";
			savename = "quadratic_0125mm";
			dx = 0.125e-3;
			order = 2;
		case 6
			type = "B";
			savename = "cubic_0125mm";
			dx = 0.125e-3;
			order = 3;
		case 7
			type = "Triangle";
			savename = "triangleLinear_0125mm";
			dx = 0.125e-3;
			order = 1;
		case 8
			type = "Triangle";
			savename = "triangleQuadratic_0125mm";
			dx = 0.125e-3;
			order = 2;
	end

	cx = 8;
	cy = 5;
	
	LX = 76.2e-3;
	LY = 25.4e-3;
	
	plotlim = 9999;
	layers = 0;
	Nx = ceil(LX/dx);
	Ny = ceil(LY/dx);
	
	surfaceMap = [0, LX; 
			  	LY, LY];
	
	meshDamageable = [12e-3 LX-12e-3;
					  0 LY-2.5e-3];

	if type=="B"
    	mesh = create_anchors_simple(order, 1, Nx, Ny, surfaceMap);
    	mesh = process_anchorsv3(mesh);

		ig = 5;
	end
	if (type=="Triangle")
		mesh = TriangleMesh(dx, LX, LY, order);

		ig = 1;
	end
	mesh.Type = type;



	mesh.nodegroups{6}.name="NoDamage";
	mesh.nodegroups{6}.nodes = [];
	for n=1:length(mesh.nodegroups{ig}.nodes)
		node = mesh.nodegroups{ig}.nodes(n);
		if (type=="B")
			coords = mesh.anchors{node}.cp;
		else
			coords = mesh.nodes(node,:);
		end
		if (coords(1)<12e-3 || coords(1)>76.2e-3-12e-3 || coords(2)>22.5e-3)
			mesh.nodegroups{6}.nodes(end+1)=node;
		end
	end


	[Nodes, NodeGroup, ElementGroup, coreData] = partition(mesh, cx, cy);
	saveToHDF(savename, Nodes, NodeGroup, ElementGroup, coreData, type)

end


function saveToHDF(fname, Nodes, NodeGroup, ElementGroup, coreData, type)
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

	clear names types;

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
		ELS = [];
		for el=1:length(ElementGroup{i}.elements)
			ELS(el,:) = uint64(ElementGroup{i}.elements{el}.cps);
		end
		h5create(fname+".h5",'/elementgroups/'+ElementGroup{i}.name,size(ELS),'Datatype','uint64');
		h5write(fname+".h5",'/elementgroups/'+ElementGroup{i}.name,ELS-1);

		if (type == "B")
			DATA = [];
			for el=1:length(ElementGroup{i}.elements)
				DATA(el,:,:) = ElementGroup{i}.elements{el}.Bezier;
			end
	
			h5create(fname+".h5",'/elementgroups/'+ElementGroup{i}.name+"_Data",size(DATA),'Datatype','double');
			h5write(fname+".h5",'/elementgroups/'+ElementGroup{i}.name+"_Data",DATA);
		end
	end

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

	h5disp(fname+".h5")
end
