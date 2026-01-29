function mesh = partition(mesh, nCores)

	coreData.NCores = nCores;

	NNodes = length(mesh.Nodes);

	if false
		NPerElem = mesh.elementgroups{1}.NNodesPerElem*mesh.elementgroups{1}.NNodesPerElem;
		ix = zeros(NPerElem*mesh.elementgroups{1}.NElems,1);
		iy = zeros(NPerElem*mesh.elementgroups{1}.NElems,1);
		for el=1:length(mesh.elementgroups{1}.elements)
			nds = mesh.elementgroups{1}.elements(el,:);
			[allocX, allocY] = ndgrid(nds,nds);
			ix(NPerElem*(el-1)+1:NPerElem*el) = allocX(:);	
			iy(NPerElem*(el-1)+1:NPerElem*el) = allocY(:);	
		end
		Connections = sparse(ix, iy, 0*ix+1, NNodes, NNodes);
	
		map = specdice(Connections,log2(nCores));
	
		clear ix iy Connections
		coreData.NodeLocs = map'+1;
	else

		xClose = [];
		yClose = [];
		zClose = [];
		for i=1:nCores/4
			xClose = [xClose mesh.params.R mesh.params.R -mesh.params.R -mesh.params.R];
			yClose = [yClose -mesh.params.R mesh.params.R -mesh.params.R mesh.params.R];
			dz = (i-0.5)*mesh.params.H/14;
			zClose = [zClose dz dz dz dz];
		end
		n = 1:length(zClose);

		F = scatteredInterpolant(xClose',yClose', zClose', n', 'nearest', 'nearest');
		coreData.NodeLocs = F(mesh.Nodes(:,1), mesh.Nodes(:,2), mesh.Nodes(:,3));
		

	end



	[~,reorderNodes] = sort(coreData.NodeLocs);

	%NodeGroups
	for eg=1:length(mesh.nodegroups)
		for el=1:length(mesh.nodegroups{eg}.nodes)
			coreData.NodeGroup{eg}.core(el) = coreData.NodeLocs(mesh.nodegroups{eg}.nodes(el));
		end
	end

	%Elements
	for eg=1:length(mesh.elementgroups)
		for el=1:size(mesh.elementgroups{eg}.elements,1)
			eNodes = mesh.elementgroups{eg}.elements(el,:);
			cNodes = coreData.NodeLocs(eNodes);
			loc = mode(cNodes);
			coreData.ElementGroup{eg}.core(el) = loc;
		end
	end

	% private/shared Nodes
	for c=1:nCores
		coreData.Nodes.PrivateNodes{c} = find(coreData.NodeLocs == c);

		coreData.Nodes.Ghosts_ImReceiving{c} = [];
		coreData.Nodes.Ghosts_ImSending{c} = [];
	end
	
	for eg=1:length(mesh.elementgroups)
		for el=1:size(mesh.elementgroups{eg}.elements,1)
			Nds = mesh.elementgroups{eg}.elements(el,:);
			NodeCores = coreData.NodeLocs(Nds);
			uniqCores = sort(unique(NodeCores));
			if (length(uniqCores)>1)
				ElemOwner = coreData.ElementGroup{eg}.core(el);
				Ghosts = find(NodeCores~=ElemOwner);
				GhostNodes = Nds(Ghosts);
				GhostCores = NodeCores(Ghosts);
				for j=1:length(Ghosts)
					coreData.Nodes.Ghosts_ImReceiving{ElemOwner}(end+1) = GhostNodes(j);
					coreData.Nodes.Ghosts_ImSending{GhostCores(j)}(end+1) = GhostNodes(j);
				end
			end
		end
	end
	for c=1:nCores
		coreData.Nodes.Ghosts_ImReceiving{c} = sort(unique(coreData.Nodes.Ghosts_ImReceiving{c}));
		coreData.Nodes.Ghosts_ImSending{c} = sort(unique(coreData.Nodes.Ghosts_ImSending{c}));
	end

	% renumbering
	mesh.Nodes = mesh.Nodes(reorderNodes,:);
	mesh.Weights = mesh.Weights(reorderNodes);
	coreData.NodeLocs = coreData.NodeLocs(reorderNodes);

	[~,OldToNew] = sort(reorderNodes);


	for ng=1:length(mesh.nodegroups)
		[~,reorderNodeGroup] = sort(coreData.NodeGroup{ng}.core);
		coreData.NodeGroup{ng}.core = coreData.NodeGroup{ng}.core(reorderNodeGroup);
		NOLD = mesh.nodegroups{ng}.nodes(reorderNodeGroup);
		mesh.nodegroups{ng}.nodes = OldToNew(NOLD);
	end

	% for ng=1:length(mesh.nodegroups)
	% 	[~,reorderNodeGroup{ng}] = sort(coreData.NodeGroup{ng}.core);
	% 	coreData.NodeGroup{ng}.core = coreData.NodeGroup{ng}.core(reorderNodeGroup{ng});
	% 	NOLD = mesh.nodegroups{ng}.nodes(reorderNodeGroup{ng});
	% 	for n=1:length(NOLD)
	% 		ToRenumber = NOLD(n);
	% 		Renumbered = find(ToRenumber==reorderNodes);
	% 		mesh.nodegroups{ng}.nodes(n) = Renumbered;
	% 		%Renumbered2 = NewIndex(NOLD(n));
	% 		%mesh.nodegroups{ng}.nodes(n) = NewIndex(NOLD(n));
	% 	end
	% end

	for eg=1:length(mesh.elementgroups)
		[~,reorderElems] = sort(coreData.ElementGroup{eg}.core);
		coreData.ElementGroup{eg}.core = coreData.ElementGroup{eg}.core(reorderElems);
		EOLD = mesh.elementgroups{eg}.elements;
		BEOLD = mesh.elementgroups{eg}.BE;

		ENEW = 0*EOLD;
		BE_NEW = 0*BEOLD;
		for n=1:size(EOLD,1)

			% for nn=1:size(EOLD,2)
			% 	ToRenumber = EOLD(reorderElems(n),nn);
			% 	Renumbered = find(ToRenumber==reorderNodes);
			% 	ENEW(n,nn) = Renumbered;
			% end

			ENEW(n,:) = OldToNew(EOLD(reorderElems(n),:));
			BE_NEW(n,:,:) = BEOLD(reorderElems(n),:,:);
		end
		mesh.elementgroups{eg}.elements = ENEW;
		mesh.elementgroups{eg}.BE = BE_NEW;
	end

	for c=1:nCores
		for n=1:length(coreData.Nodes.PrivateNodes{c})
			coreData.Nodes.PrivateNodes{c}(n) = OldToNew(coreData.Nodes.PrivateNodes{c}(n));
		end
		for n=1:length(coreData.Nodes.Ghosts_ImReceiving{c})
			coreData.Nodes.Ghosts_ImReceiving{c}(n) = OldToNew(coreData.Nodes.Ghosts_ImReceiving{c}(n));
		end
		for n=1:length(coreData.Nodes.Ghosts_ImSending{c})
			coreData.Nodes.Ghosts_ImSending{c}(n) = OldToNew(coreData.Nodes.Ghosts_ImSending{c}(n));
		end
		coreData.Nodes.PrivateNodes{c} = sort(coreData.Nodes.PrivateNodes{c});
		coreData.Nodes.Ghosts_ImReceiving{c} = sort(coreData.Nodes.Ghosts_ImReceiving{c});
		coreData.Nodes.Ghosts_ImSending{c} = sort(coreData.Nodes.Ghosts_ImSending{c});
	end
	
	for c=1:nCores
		NodeRange = find(coreData.NodeLocs==c);
		coreData.ToSave{c}.Noderange = uint64([min(NodeRange) max(NodeRange)]);
		coreData.ToSave{c}.Ghosts = uint64(coreData.Nodes.Ghosts_ImReceiving{c});

		for ng=1:length(coreData.NodeGroup)
			myNodes = find(coreData.NodeGroup{ng}.core==c);
			if isempty(myNodes)
				coreData.ToSave{c}.hasnodegroup(ng) = false;
				coreData.ToSave{c}.NodeGrouprange(ng,:) = [nan nan];
			else
				coreData.ToSave{c}.hasnodegroup(ng) = true;
				coreData.ToSave{c}.NodeGrouprange(ng,:) = uint64([min(myNodes) max(myNodes)]);
			end
		end

		for eg=1:length(coreData.ElementGroup)
			myElems = find(coreData.ElementGroup{eg}.core==c);
			if isempty(myElems)
				coreData.ToSave{c}.haselems(eg) = false;
				coreData.ToSave{c}.Elemrange(eg,:) = [nan nan];
			else
				coreData.ToSave{c}.haselems(eg) = true;
				coreData.ToSave{c}.Elemrange(eg,:) = uint64([min(myElems) max(myElems)]);
			end
		end

	end


	mesh.coreData = coreData;
end
