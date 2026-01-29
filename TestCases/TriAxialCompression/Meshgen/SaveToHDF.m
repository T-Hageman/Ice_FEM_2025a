function SaveToHDF(fname, mesh)
	if exist(fname+".h5", 'file')==2
		delete(fname+".h5");
	end

	% nodes
	h5create(fname+".h5",'/nodes',size(mesh.Nodes),'Datatype','double');
	h5write(fname+".h5",'/nodes', mesh.Nodes)

	% NodeData
	h5create(fname+".h5",'/node_DataTypes',1,'Datatype','string');
	h5write(fname+".h5",'/node_DataTypes', "W")

	h5create(fname+".h5",'/node_Data',size(mesh.Weights),'Datatype','double');
	h5write(fname+".h5",'/node_Data', mesh.Weights)

	% nodegroups
	h5create(fname+".h5",'/nodegroupnames',length(mesh.nodegroups),'Datatype','string');
	for i=1:length(mesh.nodegroups)
		names{i} = mesh.nodegroups{i}.name;
	end

	h5write(fname+".h5",'/nodegroupnames', string(names))
	for i=1:length(mesh.nodegroups)
		h5create(fname+".h5",'/nodegroups/'+mesh.nodegroups{i}.name,length(mesh.nodegroups{i}.nodes),'Datatype','uint64');
		h5write(fname+".h5",'/nodegroups/'+mesh.nodegroups{i}.name,uint64(mesh.nodegroups{i}.nodes)-1);
	end

	% Elementgroups
	clear names
	h5create(fname+".h5",'/elementgroupnames',length(mesh.elementgroups),'Datatype','string');
	h5create(fname+".h5",'/elementgrouptypes',length(mesh.elementgroups),'Datatype','string');
	for i=1:length(mesh.elementgroups)
		names{i} = mesh.elementgroups{i}.name;
		types{i} = mesh.elementgroups{i}.type;
	end

	h5write(fname+".h5",'/elementgroupnames', string(names))
	h5write(fname+".h5",'/elementgrouptypes', string(types))
	for i=1:length(mesh.elementgroups)
		ELS = uint64(mesh.elementgroups{i}.elements-1); %correct for starting at zero
		DATA = mesh.elementgroups{i}.BE;

		h5create(fname+".h5",'/elementgroups/'+mesh.elementgroups{i}.name,size(ELS),'Datatype','uint64');
		h5write(fname+".h5",'/elementgroups/'+mesh.elementgroups{i}.name,ELS);

		h5create(fname+".h5",'/elementgroups/'+mesh.elementgroups{i}.name+"_Data",size(DATA),'Datatype','double');
		h5write(fname+".h5",'/elementgroups/'+mesh.elementgroups{i}.name+"_Data",DATA);
	end

	% partition data
	h5create(fname+".h5",'/NCores', 1,'Datatype','uint64');
	h5write(fname+".h5",'/NCores', uint64(mesh.coreData.NCores));

	for c=1:mesh.coreData.NCores
		h5create(fname+".h5",'/Partition/'+string(c-1)+'/noderange', 2,'Datatype','uint64');
		h5write(fname+".h5",'/Partition/'+string(c-1)+'/noderange',uint64(mesh.coreData.ToSave{c}.Noderange-1));

		if (isempty(mesh.coreData.ToSave{c}.Ghosts))
			h5create(fname+".h5",'/Partition/'+string(c-1)+'/Hasghosts', 1,'Datatype','uint8');
			h5write(fname+".h5",'/Partition/'+string(c-1)+'/Hasghosts', uint8(false));
		else
			h5create(fname+".h5",'/Partition/'+string(c-1)+'/Hasghosts', 1,'Datatype','uint8');
			h5write(fname+".h5",'/Partition/'+string(c-1)+'/Hasghosts', uint8(true));

			h5create(fname+".h5",'/Partition/'+string(c-1)+'/ghosts', length(mesh.coreData.ToSave{c}.Ghosts),'Datatype','uint64');
			h5write(fname+".h5",'/Partition/'+string(c-1)+'/ghosts',uint64(mesh.coreData.ToSave{c}.Ghosts-1));
		end

		h5create(fname+".h5",'/Partition/'+string(c-1)+'/Haselems', length(mesh.coreData.ToSave{c}.haselems),'Datatype','uint8');
		h5write(fname+".h5",'/Partition/'+string(c-1)+'/Haselems',uint8(mesh.coreData.ToSave{c}.haselems));

		h5create(fname+".h5",'/Partition/'+string(c-1)+'/elemrange', size(mesh.coreData.ToSave{c}.Elemrange),'Datatype','uint64');
		h5write(fname+".h5",'/Partition/'+string(c-1)+'/elemrange',uint64(mesh.coreData.ToSave{c}.Elemrange-1));

		h5create(fname+".h5",'/Partition/'+string(c-1)+'/hasnodegroup', length(mesh.coreData.ToSave{c}.hasnodegroup),'Datatype','uint8');
		h5write(fname+".h5",'/Partition/'+string(c-1)+'/hasnodegroup',uint8(mesh.coreData.ToSave{c}.hasnodegroup));

		h5create(fname+".h5",'/Partition/'+string(c-1)+'/nodegrouprange', size(mesh.coreData.ToSave{c}.NodeGrouprange),'Datatype','uint64');
		h5write(fname+".h5",'/Partition/'+string(c-1)+'/nodegrouprange',uint64(mesh.coreData.ToSave{c}.NodeGrouprange-1));
	end

	%h5disp(fname+".h5")
end