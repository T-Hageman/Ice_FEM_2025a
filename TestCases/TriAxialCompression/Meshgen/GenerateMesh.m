close all
clear
clc

addpath(genpath('Z:\Libraries\igafem'))
addpath(genpath('Z:\Libraries\meshpart-master'))

tic

ncores = 56;
for dx=[10, 5, 2.5, 2, 1]
	for p=[2]
		% Mesh parameters
		mesh.params.p = p;
		mesh.params.C0 = false;
		mesh.params.R = 10e-2/2;
		mesh.params.H = 30e-2;
		mesh.params.dx = dx*1e-3;
		
		filename = "Mesh_"+string(mesh.params.dx*1000)+"mm_p"+string(mesh.params.p)+"_c"+string(ncores);
		fprintf(filename+"\n");
		
		% element sizes
		mesh.params.Nr = ceil(mesh.params.R/mesh.params.dx);
		mesh.params.Ny = ceil(mesh.params.H/mesh.params.dx);
		
		% basic geometry
		fprintf("\tCreating geometry\n")
		mesh = GenerateCylinder(mesh);
		fprintf("\tfinished\n")
		
		% mesh
		fprintf("\tCreating mesh\n")
		mesh = CreateMesh(mesh);
		fprintf("\tfinished\n")
		
		% Partition
		fprintf("\tPartitioning mesh\n")
		mesh = partition(mesh, ncores);
		fprintf("\tfinished\n")
		
		% save to file
		fprintf("\tSaving mesh\n")
		SaveToHDF(filename, mesh)
		fprintf("\tfinished\n")
		
		% plotting
		fprintf("\tPlotting\n")
		PrintInfo(mesh);
		%PlotMesh(mesh, true);
		PlotPartition(mesh);
		fprintf("\tfinished\n")
		
		
		toc
		clear mesh
		close all
	end
end



function PrintInfo(mesh)
	fprintf("Mesh Information\n");

	fprintf("\t Total number of Nodes: "+string(size(mesh.Nodes,1))+"\n");
	for i=1:length(mesh.nodegroups)
		fprintf("\t\t "+mesh.nodegroups{i}.name+": ");
		fprintf(string(length(mesh.nodegroups{i}.nodes))+"\n");
	end

	fprintf("\t Number of Elements: "+"\n");
	for i=1:length(mesh.elementgroups)
		fprintf("\t\t "+mesh.elementgroups{i}.name+": ");
		fprintf(string(size(mesh.elementgroups{i}.elements,1))+"  ("+mesh.elementgroups{i}.type+")\n");
	end

	fprintf("\t Other statistics:\n")
	fprintf("\t\t Volume: %e \n", (pi*mesh.params.R^2*mesh.params.H))
	fprintf("\t\t Bottom Area: %e \n", (pi*mesh.params.R^2))
	fprintf("\t\t Side Area: %e \n", (2*pi*mesh.params.R*mesh.params.H))
end

function PlotMesh(mesh, PlotNodes)
	c = {"b","c","m","r","k","b","c","m"};

	figure
	nrbkntplot(mesh.Nb);
	hold on
	viscircles([0,0],mesh.params.R);

	%nodes
	if (PlotNodes)
		%plot3(mesh.Nodes(:,1),mesh.Nodes(:,2),mesh.Nodes(:,3),'r*')
		for i=1:length(mesh.nodegroups)
			plot3(mesh.Nodes(mesh.nodegroups{i}.nodes,1),mesh.Nodes(mesh.nodegroups{i}.nodes,2),mesh.Nodes(mesh.nodegroups{i}.nodes,3),c{i}+'*')
		end
	end

	axis equal
	view(3);
end

function PlotPartition(mesh)
	clr = distinguishable_colors(mesh.coreData.NCores);

	figure
	for c=1:mesh.coreData.NCores
		N = mesh.Nodes(mesh.coreData.NodeLocs==c,:);
		plot3(N(:,1), N(:,2), N(:,3), '*','Color',clr(c,:));
		hold on
	end

	axis equal
	view(3);
end

function mesh = GenerateCylinder(mesh)
	dim = 3;
	q = 3; %Number of points in U
	r = 3; %Number of points in V
	s = 3; %Number of points in W
	
	uKnot = [0 0 0 1 1 1];
	vKnot = [0 0 0 1 1 1];
	wKnot = [0 0 0 1 1 1];
	knots  = {uKnot vKnot wKnot};
	
	sinR = sind(45)*mesh.params.R;
	dC = mesh.params.R*sqrt(2);
	W = sqrt(2)/2;
	
	ControlP_in = zeros(dim+1,q,r,s); %for constructing NURBS
	
	ControlP_in(dim+1,:,:,:)   = 1;
	for i=1:3
		Z = mesh.params.H*(i-1)/2;
		ControlP_in(1:3,1,1,i) = [-sinR;-sinR;Z];
		ControlP_in(1:3,2,1,i) = [0;-dC;Z];
		ControlP_in(1:3,3,1,i) = [sinR;-sinR;Z];
		ControlP_in(1:3,1,2,i) = [-dC;0;Z];
		ControlP_in(1:3,2,2,i) = [0;0;Z];
		ControlP_in(1:3,3,2,i) = [dC;0;Z];
		ControlP_in(1:3,1,3,i) = [-sinR;sinR;Z];
		ControlP_in(1:3,2,3,i) = [0;dC;Z];
		ControlP_in(1:3,3,3,i) = [sinR;sinR;Z];
	
		ControlP_in(dim+1,2,1,i)   = W;
		ControlP_in(dim+1,1,2,i)   = W;
		ControlP_in(dim+1,3,2,i)   = W;
		ControlP_in(dim+1,2,3,i)   = W;
	end
	
	ControlP_in(1,:,:,:) = ControlP_in(1,:,:,:).*ControlP_in(dim+1,:,:,:);
	ControlP_in(2,:,:,:) = ControlP_in(2,:,:,:).*ControlP_in(dim+1,:,:,:);
	ControlP_in(3,:,:,:) = ControlP_in(3,:,:,:).*ControlP_in(dim+1,:,:,:);
	mesh.Nb = nrbmak(ControlP_in,knots);
end

function mesh = CreateMesh(mesh)
	%mesh order & insert elements
	Order_Increase = mesh.params.p-2;
	if (Order_Increase>0)
		mesh.Nb = nrbdegelev(mesh.Nb, [Order_Increase,Order_Increase,Order_Increase]);
	end

	xi_Kinsert = zeros(1,2*mesh.params.Nr-1);
	for i=1:2*mesh.params.Nr-1
		xi_Kinsert(i)=i/(2*mesh.params.Nr);
	end
	
	eta_Kinsert = zeros(1,mesh.params.Ny-1);
	for i=1:mesh.params.Ny-1
		eta_Kinsert(i)=i/(mesh.params.Ny);
	end  

	mesh.Nb = nrbkntins(mesh.Nb,{xi_Kinsert xi_Kinsert eta_Kinsert});
	if (mesh.params.C0==true)
		mesh.Nb = nrbkntins(mesh.Nb,{xi_Kinsert xi_Kinsert eta_Kinsert});
	end

	% Control points
	uo     = mesh.Nb.order(1)-1;
    vo     = mesh.Nb.order(2)-1;
	wo	   = mesh.Nb.order(3)-1;
    uKnot  = cell2mat(mesh.Nb.knots(1));
    vKnot  = cell2mat(mesh.Nb.knots(2));
	wKnot  = cell2mat(mesh.Nb.knots(3));
    noPtsX = length(uKnot)-uo-1;
    noPtsY = length(vKnot)-vo-1;
	noPtsZ = length(wKnot)-wo-1;

	Nds(:,:,:,:)= mesh.Nb.coefs(:,:,:,:);

	count = 0;
	mesh.Nodes = zeros(noPtsX*noPtsY*noPtsZ,3);
	mesh.Weights = zeros(noPtsX*noPtsY*noPtsZ,1);
	NodeNums = zeros(noPtsX, noPtsY, noPtsZ);
    for k=1:noPtsZ
        for j=1:noPtsY
			for i=1:noPtsX
				count = count + 1;
            	NodeNums(i,j,k) = count;
				mesh.Nodes(count,:) = Nds(1:3,i,j,k)/Nds(4,i,j,k);
				mesh.Weights(count) = Nds(4,i,j,k);
			end
        end
    end

	% Node Groups
	mesh.nodegroups{1}.name = "bottom";
	mesh.nodegroups{1}.nodes= reshape(NodeNums(:,:,1),1,[])';

    mesh.nodegroups{2}.name = "top";
	mesh.nodegroups{2}.nodes= reshape(NodeNums(:,:,end),1,[])';

    mesh.nodegroups{3}.name = "side";
	SideNodes = [reshape(NodeNums(1,:,:),1,[])';
				reshape(NodeNums(end,:,:),1,[])';
				reshape(NodeNums(:,1,:),1,[])';
				reshape(NodeNums(:,end,:),1,[])'];

	mesh.nodegroups{3}.nodes = SideNodes;

	mesh.nodegroups{4}.name = "bottom_HorizontalLine";
	mesh.nodegroups{4}.nodes=[];
	mesh.nodegroups{5}.name = "bottom_VerticalLine";
	mesh.nodegroups{5}.nodes=[];

	dx_max = mesh.params.dx;
	for i=1:length(mesh.nodegroups{1}.nodes)
		coords = mesh.Nodes(mesh.nodegroups{1}.nodes(i),:);

		if (abs(coords(1))<dx_max)
			mesh.nodegroups{5}.nodes=[mesh.nodegroups{5}.nodes; mesh.nodegroups{1}.nodes(i)];
		end
		if (abs(coords(2))<dx_max)
			mesh.nodegroups{4}.nodes=[mesh.nodegroups{4}.nodes; mesh.nodegroups{1}.nodes(i)];
		end
	end

	mesh.nodegroups{6}.name = "top_HorizontalLine";
	mesh.nodegroups{6}.nodes=[];
	mesh.nodegroups{7}.name = "top_VerticalLine";
	mesh.nodegroups{7}.nodes=[];

	dx_max = mesh.params.dx;
	for i=1:length(mesh.nodegroups{2}.nodes)
		coords = mesh.Nodes(mesh.nodegroups{2}.nodes(i),:);

		if (abs(coords(1))<dx_max)
			mesh.nodegroups{7}.nodes=[mesh.nodegroups{7}.nodes; mesh.nodegroups{2}.nodes(i)];
		end
		if (abs(coords(2))<dx_max)
			mesh.nodegroups{6}.nodes=[mesh.nodegroups{6}.nodes; mesh.nodegroups{2}.nodes(i)];
		end
	end

	mesh.nodegroups{8}.name = "internal";
	mesh.nodegroups{8}.nodes= 1:noPtsX*noPtsY*noPtsZ;

	% Element Groups and Bezier extractors
	uKnotVec = unique(uKnot);
    vKnotVec = unique(vKnot);
	wKnotVec = unique(wKnot);

    noElemsU = length(uKnotVec)-1; % # of elements xi dir.
    noElemsV = length(vKnotVec)-1; % # of elements eta dir.
	noElemsW = length(wKnotVec)-1; % # of elements eta dir.

	[~,elConnU] = buildConnectivity(uo,uKnot,noElemsU);
    [~,elConnV] = buildConnectivity(vo,vKnot,noElemsV);
	[~,elConnW] = buildConnectivity(wo,wKnot,noElemsW);

	%Interior
	mesh.elementgroups{1}.name = "internal";
	mesh.elementgroups{1}.type = "NURBS3_Cube"+string(uo);
	mesh.elementgroups{1}.NElems = noElemsU*noElemsV*noElemsW;
	mesh.elementgroups{1}.NNodesPerElem = size(elConnU,2)*size(elConnV,2)*size(elConnW,2);

	C = bezierExtraction3D(uKnot, vKnot, wKnot,uo,vo,wo);

	Elems = zeros(mesh.elementgroups{1}.NElems, mesh.elementgroups{1}.NNodesPerElem);
	BE = zeros(mesh.elementgroups{1}.NElems, mesh.elementgroups{1}.NNodesPerElem, mesh.elementgroups{1}.NNodesPerElem);
	counter = 0;
	for w=1:noElemsW
		for v=1:noElemsV
			for u=1:noElemsU
				counter = counter + 1;
				Elems(counter, :) = reshape(NodeNums(elConnU(u,:),elConnV(v,:),elConnW(w,:)),1,[]);
				BE(counter,:,:) = C{u,v,w}(:,:);
			end
		end
	end
	mesh.elementgroups{1}.elements = Elems;
	mesh.elementgroups{1}.BE = BE;

	%bottom/top
	mesh.elementgroups{2}.name = "bottom";
	mesh.elementgroups{3}.name = "top";
	mesh.elementgroups{2}.type = "NURBS3_Plane"+string(uo);
	mesh.elementgroups{3}.type = "NURBS3_Plane"+string(uo);

	mesh.elementgroups{2}.NElems = noElemsU*noElemsV;
	mesh.elementgroups{3}.NElems = noElemsU*noElemsV;

	mesh.elementgroups{2}.NNodesPerElem = size(elConnU,2)*size(elConnV,2);
	mesh.elementgroups{3}.NNodesPerElem = size(elConnU,2)*size(elConnV,2);

	C = bezierExtraction2D(uKnot, vKnot,uo,vo);

	ElemsBot = zeros(mesh.elementgroups{2}.NElems, mesh.elementgroups{2}.NNodesPerElem);
	ElemsTop = zeros(mesh.elementgroups{3}.NElems, mesh.elementgroups{3}.NNodesPerElem);
	BEBot = zeros(mesh.elementgroups{2}.NElems, mesh.elementgroups{2}.NNodesPerElem, mesh.elementgroups{2}.NNodesPerElem);
	BETop = zeros(mesh.elementgroups{3}.NElems, mesh.elementgroups{3}.NNodesPerElem, mesh.elementgroups{3}.NNodesPerElem);
	counter = 0;
	for v=1:noElemsV
		for u=1:noElemsU
			counter = counter + 1;
			ElemsBot(counter, :) = reshape(NodeNums(elConnU(u,:),elConnV(v,:),1),1,[]);
			ElemsTop(counter, :) = reshape(NodeNums(elConnU(u,:),elConnV(v,:),end),1,[]);
			BEBot(counter,:,:) = C{u,v}(:,:);
			BETop(counter,:,:) = C{u,v}(:,:);
		end
	end
	mesh.elementgroups{2}.elements = ElemsBot;
	mesh.elementgroups{2}.BE = BEBot;	
	mesh.elementgroups{3}.elements = ElemsTop;
	mesh.elementgroups{3}.BE = BETop;
	
	%side
    mesh.elementgroups{4}.name = "side";
	mesh.elementgroups{4}.type = "NURBS3_Plane"+string(uo);
	mesh.elementgroups{4}.NElems = 2*noElemsU*noElemsW+2*noElemsV*noElemsW;
	mesh.elementgroups{4}.NNodesPerElem = size(elConnU,2)*size(elConnW,2);

	C = bezierExtraction2D(uKnot, wKnot,uo,wo);
	ElemsFront = zeros(noElemsU*noElemsW, mesh.elementgroups{4}.NNodesPerElem);
	ElemsBack = zeros(noElemsU*noElemsW, mesh.elementgroups{4}.NNodesPerElem);
	BEFront = zeros(noElemsU*noElemsW, mesh.elementgroups{4}.NNodesPerElem, mesh.elementgroups{4}.NNodesPerElem);
	BEBack = zeros(noElemsU*noElemsW, mesh.elementgroups{4}.NNodesPerElem, mesh.elementgroups{4}.NNodesPerElem);
	counter = 0;
	for w=1:noElemsW
		for u=1:noElemsU
			counter = counter + 1;
			ElemsFront(counter, :) = reshape(NodeNums(elConnU(u,:),1,elConnW(w,:)),1,[]);
			ElemsBack(counter, :) = reshape(NodeNums(elConnU(u,:),end,elConnW(w,:)),1,[]);
			BEFront(counter,:,:) = C{u,w}(:,:);
			BEBack(counter,:,:) = C{u,w}(:,:);
		end
	end

	C = bezierExtraction2D(vKnot, wKnot,vo,wo);
	ElemsLeft = zeros(noElemsV*noElemsW, mesh.elementgroups{4}.NNodesPerElem);
	Elemsright = zeros(noElemsV*noElemsW, mesh.elementgroups{4}.NNodesPerElem);
	BELeft = zeros(noElemsV*noElemsW, mesh.elementgroups{4}.NNodesPerElem, mesh.elementgroups{4}.NNodesPerElem);
	BERight = zeros(noElemsV*noElemsW, mesh.elementgroups{4}.NNodesPerElem, mesh.elementgroups{4}.NNodesPerElem);
	counter = 0;
	for w=1:noElemsW
		for v=1:noElemsV
			counter = counter + 1;
			ElemsLeft(counter, :) = reshape(NodeNums(1,elConnV(v,:),elConnW(w,:)),1,[]);
			Elemsright(counter, :) = reshape(NodeNums(end,elConnV(v,:),elConnW(w,:)),1,[]);
			BELeft(counter,:,:) = C{v,w}(:,:);
			BERight(counter,:,:) = C{v,w}(:,:);
		end
	end
	mesh.elementgroups{4}.elements = [ElemsFront; ElemsBack; ElemsLeft; Elemsright];
	mesh.elementgroups{4}.BE = [BEFront; BEBack; BELeft; BERight];
end


