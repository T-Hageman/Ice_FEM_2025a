function [Nodes, NodeGroup, ElementGroup, coreData] = partition(mesh, cx, cy)
	Nodes = [];
	for i=1:length(mesh.anchors)
		Nodes(i,:) = mesh.anchors{i}.cp;
	end
	coreData.NCores = cx*cy;

	%% allocation
	dx = (max(Nodes(:,1))-min(Nodes(:,1)))/cx; dxmn = min(Nodes(:,1));
	dy = max(Nodes(:,2))/cy;

	NearPoints = 1:cx*cy;
	NearPoints_Loc = [];
	for i=1:cx
		for j=1:cy
			point(1) = (i-0.5)*dx+1e-6*rand+dxmn;
			point(2) = (j-0.5)*dy+1e-6*rand;
			NearPoints_Loc(end+1,1:2) = point;
		end
	end

	%Nodes
	Allocator = scatteredInterpolant(NearPoints_Loc(:,1),NearPoints_Loc(:,2),NearPoints','nearest');
	coreData.NodeLocs = Allocator(Nodes(:,1),Nodes(:,2));

	%NodeGroups
	NodeGroup = mesh.nodegroups;
	for eg=1:length(NodeGroup)
		for el=1:length(NodeGroup{eg}.nodes)
			coreData.NodeGroup{eg}.core(el) = coreData.NodeLocs(NodeGroup{eg}.nodes(el));
		end
	end

	%Elements
	ElementGroup = mesh.elementgroups;
	for eg=1:length(ElementGroup)
		for el=1:length(ElementGroup{eg}.elements)
			elMean = mean(Nodes(ElementGroup{eg}.elements{el}.cps,:));
			coreData.ElementGroup{eg}.core(el) = Allocator(elMean);
		end
	end

	% private/shared Nodes
	for c=1:cx*cy
		coreData.Nodes.PrivateNodes{c} = find(coreData.NodeLocs == c);

		coreData.Nodes.Ghosts_ImReceiving{c} = [];
		coreData.Nodes.Ghosts_ImSending{c} = [];
	end
	
	for eg=1:length(ElementGroup)
		for el=1:length(ElementGroup{eg}.elements)
			Nds = ElementGroup{eg}.elements{el}.cps;
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
	for c=1:cx*cy
		coreData.Nodes.Ghosts_ImReceiving{c} = sort(unique(coreData.Nodes.Ghosts_ImReceiving{c}));
		coreData.Nodes.Ghosts_ImSending{c} = sort(unique(coreData.Nodes.Ghosts_ImSending{c}));
	end

	% renumbering
	[~,reorderNodes] = sort(coreData.NodeLocs);
	Nodes = Nodes(reorderNodes,:);
	coreData.NodeLocs = coreData.NodeLocs(reorderNodes);

	for ng=1:length(NodeGroup)
		[~,reorderNodeGroup{ng}] = sort(coreData.NodeGroup{ng}.core);
		coreData.NodeGroup{ng}.core = coreData.NodeGroup{ng}.core(reorderNodeGroup{ng});
		NOLD = NodeGroup{ng}.nodes;
		for n=1:length(NOLD)
			ToRenumber = NOLD(reorderNodeGroup{ng}(n));
			Renumbered = find(ToRenumber==reorderNodes);
			NodeGroup{ng}.nodes(n) = Renumbered;
		end
	end

	for eg=1:length(ElementGroup)
		[~,reorderElems{eg}] = sort(coreData.ElementGroup{eg}.core);
		coreData.ElementGroup{eg}.core = coreData.ElementGroup{eg}.core(reorderElems{eg});
		EOLD = ElementGroup{eg}.elements;
		
		for n=1:length(EOLD)
			ElementGroup{eg}.elements{n} = EOLD{reorderElems{eg}(n)};
			for nn=1:length(EOLD{1}.cps)
				ToRenumber = EOLD{reorderElems{eg}(n)}.cps(nn);
				Renumbered = find(ToRenumber==reorderNodes);
				ElementGroup{eg}.elements{n}.cps(nn) = Renumbered;
			end
		end
	end

	for c=1:cx*cy
		for n=1:length(coreData.Nodes.PrivateNodes{c})
			coreData.Nodes.PrivateNodes{c}(n) = find(coreData.Nodes.PrivateNodes{c}(n)==reorderNodes);
		end
		for n=1:length(coreData.Nodes.Ghosts_ImReceiving{c})
			coreData.Nodes.Ghosts_ImReceiving{c}(n) = find(coreData.Nodes.Ghosts_ImReceiving{c}(n)==reorderNodes);
		end
		for n=1:length(coreData.Nodes.Ghosts_ImSending{c})
			coreData.Nodes.Ghosts_ImSending{c}(n) = find(coreData.Nodes.Ghosts_ImSending{c}(n)==reorderNodes);
		end
		coreData.Nodes.PrivateNodes{c} = sort(coreData.Nodes.PrivateNodes{c});
		coreData.Nodes.Ghosts_ImReceiving{c} = sort(coreData.Nodes.Ghosts_ImReceiving{c});
		coreData.Nodes.Ghosts_ImSending{c} = sort(coreData.Nodes.Ghosts_ImSending{c});
	end
	
	for c=1:cx*cy
		NodeRange = find(coreData.NodeLocs==c);
		coreData.ToSave{c}.Noderange = uint64([min(NodeRange) max(NodeRange)]);
		coreData.ToSave{c}.Ghosts = uint64(coreData.Nodes.Ghosts_ImReceiving{c});

		for ng=1:length(NodeGroup)
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

	plotMesh(Nodes, ElementGroup, coreData, false);
end

function plotMesh(nodes, ElementGroup, coreData, plottext)
	clrs = distinguishable_colors(coreData.NCores);

	figure
	for i=1:coreData.NCores
		plot(nodes(coreData.NodeLocs==i,1),nodes(coreData.NodeLocs==i,2),'*','Color',clrs(i,:));
		hold on
	end
	if (plottext)
		text(nodes(:,1),nodes(:,2),string(1:length(nodes)));
	end
	axis equal

	% for eg=1:length(ElementGroup)/2
	% 	plotElems(ElementGroup{eg}, nodes, coreData.ElementGroup{eg}, clrs)
	% end


% 	subplot(2,1,2)
% 	nds = unique(ElementGroup{6}.Elements);
% 	nlocs = coreData.NodeLocs(nds);
% 
% 	for i=1:coreData.NCores
% 		plot(nodes(nds(nlocs==i),1),nodes(nds(nlocs==i),2),'*','Color',clrs(i,:))
% 		hold on
% 	end
% 	if (plottext)
% 		text(nodes(nds,1),nodes(nds,2),string(nds));
% 	end
% 	for eg=length(ElementGroup)/2+1:length(ElementGroup)
% 		plotElems(ElementGroup{eg}, nodes, coreData.ElementGroup{eg}, clrs)
% 	end
end

function plotElems(ElementGroup, nodes, coreData, clrs)
	if (ElementGroup.type == "B2_3" || ElementGroup.type == "B2_4")
		for el=1:length(ElementGroup.elements)
			crns = ElementGroup.elements{el}.corners_phys;
			plot([crns(1,1) crns(1,2) crns(1,2) crns(1,1) crns(1,1)]' ,[crns(2,1) crns(2,1) crns(2,2) crns(2,2) crns(2,1)]','Color',clrs(coreData.core(el),:));
		end
	end
end


