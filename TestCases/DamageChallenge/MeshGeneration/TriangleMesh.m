function mesh = TriangleMesh(dx, LX, LY, order)

	R_1 = [3,4,0,LX,LX,0,0,0,LY,LY]';
	gm = [R_1];
	sf = 'R1';
	ns = char('R1');
	ns = ns';
	shape = decsg(gm,sf,ns);
	
	geo = createpde(1);
	geometryFromEdges(geo,shape);
	generateMesh(geo,'Hgrad',1.2,'Hface',{1,dx});

	% pdegplot(geo,"VertexLabels","on","EdgeLabels","on","FaceLabels","on")
	% 
	% figure
	% tiledlayout('flow')
	% nexttile
	% pdegplot(geo,"VertexLabels","on","EdgeLabels","on","FaceLabels","on")
	% nexttile
	% pdeplot(geo,'NodeLabels','off','ElementLabels','off')

	Nodes = geo.Mesh.Nodes';
	%% interior elements
		Elementgroups{1}.Name = "internal";
		Elementgroups{1}.Type = "T6B";
		Elementgroups{1}.elements = geo.Mesh.Elements(:,findElements(geo.Mesh,'region','Face',1))';
	
	%% exterior boundary elements
		Elementgroups{2}.Name = "left";
			Elementgroups{2}.Type = "L3B";
			N = findNodes(geo.Mesh,'region','Edge',[4]);
			xy = Nodes(N,:,:);
			[~,i] = sort(xy(:,2));
			cntr = 0;
			while length(i)>2
				cntr = cntr+1;
				Elementgroups{2}.elements(cntr,:) = N(i(1:3));
				i(1:2) = [];
			end
	
		Elementgroups{3}.Name = "right";
			Elementgroups{3}.Type = "L3B";
			N = findNodes(geo.Mesh,'region','Edge',[2]);
			xy = Nodes(N,:,:);
			[~,i] = sort(xy(:,2));
			cntr = 0;
			while length(i)>2
				cntr = cntr+1;
				Elementgroups{3}.elements(cntr,:) = N(i(1:3));
				i(1:2) = [];
			end
	
		Elementgroups{4}.Name = "top";
			Elementgroups{4}.Type = "L3B";
			N = findNodes(geo.Mesh,'region','Edge',[3]);
			xy = Nodes(N,:,:);
			[~,i] = sort(xy(:,1));
			cntr = 0;
			while length(i)>2
				cntr = cntr+1;
				Elementgroups{4}.elements(cntr,:) = N(i(1:3));
				i(1:2) = [];
			end
	
		Elementgroups{5}.Name = "bottom";
			Elementgroups{5}.Type = "L3B";
			N = findNodes(geo.Mesh,'region','Edge',[1]);
			xy = Nodes(N,:,:);
			[~,i] = sort(xy(:,1));
			cntr = 0;
			while length(i)>2
				cntr = cntr+1;
				Elementgroups{5}.elements(cntr,:) = N(i(1:3));
				i(1:2) = [];
			end
	
	
	if (order == 1)
		n_groups = length(Elementgroups);
		for g=1:n_groups
			groupIdx = g;
			Elementgroups{groupIdx}.Name = Elementgroups{g}.Name;
			Elementgroups{g}.Name = Elementgroups{g}.Name;
		
			if (g==1)
				Elementgroups{groupIdx}.Type = "T3B";
				Elementgroups{groupIdx}.elements = Elementgroups{g}.elements(:,1:3);
			else
				Elementgroups{groupIdx}.Type = "L2B";
				Elementgroups{groupIdx}.elements = Elementgroups{g}.elements(:,[1,3]);
			end
		
		end
	end
	
	
	for g=1:length(Elementgroups)
    	Nodegroups{g}.name = Elementgroups{g}.Name;
    	Nodegroups{g}.nodes = unique(reshape(Elementgroups{g}.elements,[],1));
	end

	mesh.nodes = Nodes;
	for g=1:length(Elementgroups)
		mesh.elementgroups{g}.name = Elementgroups{g}.Name;
		mesh.elementgroups{g}.type = Elementgroups{g}.Type;
		for i=1:length(Elementgroups{g}.elements)
			mesh.elementgroups{g}.elements{i}.cps = Elementgroups{g}.elements(i,:);
		end
	end
	mesh.nodegroups = Nodegroups;
end