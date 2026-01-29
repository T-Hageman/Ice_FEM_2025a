%% t-splines unequal order mesh generation
close all
addpath(genpath('../../../../Libraries/igafem'))
clearvars
clc
	delete(gcp('nocreate'))
	parpool('threads')


cx = 4;
cy = 14;

for dx=[2]*1e-3
	for p=[2,3]
		LX = 10e-2;
		LY = 30.0e-2;
		
		plotlim = 9999;
		layers = 0;
		Nx = ceil(LX/dx)
		Ny = ceil(LY/dx)
		
		savename = "Plate2D_p"+string(p)+"_"+string(dx)+"m_"+string(cx*cy)+"c";
		surfaceMap = [-LX/2, LX/2; 
		  			LY, LY];
		
		if true
    		mesh1 = create_anchors_simple(p, 1, Nx, Ny, surfaceMap);
    		%mesh2 = create_anchors_simple(3, 2, Nx, Ny, surfaceMap);
		
    		mesh1 = process_anchorsv3(mesh1)
    		%mesh2 = process_anchorsv3(mesh2)
    		%save("mesh_"+savename+".mat","mesh1","mesh2")
		else
    		%load("mesh_"+savename+".mat")
		end
		%[mesh_total, n_off1, n_off2] = combine_meshes(mesh1, mesh2, "2_", "3_");
		mesh_total = mesh1;
		
		[Nodes, NodeGroup, ElementGroup, coreData] = partition(mesh_total, cx, cy);
		saveToHDF(savename, Nodes, NodeGroup, ElementGroup, coreData)
		
		
		n_off1 = 1;
		n_off2 = 1;
		
		%% plot parametric domain
		%f1 = figure('Name','Parametric domain');
		%subplot(1,2,1)
		%plot_par_domain(mesh1, n_off1-1, plotlim)
		%box off
		
		%subplot(1,2,2)
		%plot_par_domain(mesh2, n_off2-1, plotlim)
		%box off
		
		%% plot physical domain
		% f3 = figure('Name','Physical domain');
		% subplot(1,2,1)
		% plot_phys_domain(mesh1, n_off1-1)
		
		%subplot(1,2,2)
		%plot_phys_domain(mesh2, n_off2-1)
	end
end

%% related functions
function plot_par_domain(mesh, offset, plotlim)
    
    %% anchor locations
    pl=0;
    for an = 1:length(mesh.anchors)  
        anchor(1) = median(mesh.anchors{an}.knot_ind_x);
        anchor(2) = median(mesh.anchors{an}.knot_ind_y);
        if (anchor(1)<plotlim)
            pl=pl+1;
            plot_anchor(pl,1:2) = anchor;
            %text(anchor(1),anchor(2),int2str(an+offset),'FontSize',8)    
        end
        hold on     
    end
    plot(plot_anchor(:,1), plot_anchor(:,2),'k*');
    
    %% reduced continuity lines
    rp=0;
    for an=1:length(mesh.anchors)  %% reduced cont. lines
        knot_x = mesh.anchors{an}.knot_ind_x;
        knot_y = mesh.anchors{an}.knot_ind_y;
        
        for i=1:length(knot_x)-1
            for j=1:length(knot_y)-1
                corners = [knot_x(i), knot_y(j)];
                len = [knot_x(i+1)-knot_x(i), knot_y(j+1)-knot_y(j)];
                pos = [corners, len];
                
                if (len(1)~=0 && len(2)~=0)
                    if corners(1)<plotlim
                        rp=rp+1;
                        x_rec(rp,1:4) = [knot_x(i) knot_x(i+1) knot_x(i+1), knot_x(i)];
                        y_rec(rp,1:4) = [knot_y(j) knot_y(j) knot_y(j+1) knot_y(j+1)];
                        
                        rectangle('Position',pos,'EdgeColor','r')
                    end
                end
                
            end
        end 
    end
    patch('XData',x_rec', 'YData',y_rec', 'EdgeColor','red','FaceColor','none')
    
    %% meshlines
    if (mod(mesh.p,2)==0) %% even mesh order
        plot_meshlines_even(mesh, plotlim);
    else %odd mesh order
        plot_meshlines_odd(mesh, plotlim);
    end
    
    if (true)  
        %% elements
        clear x_rec y_rec
        rp = 0;
        for el=1:length(mesh.elementgroups{1}.elements)
            crn = mesh.elementgroups{1}.elements{el}.corners_ind;

            x = [crn(1,1) crn(1,2) crn(1,2) crn(1,1)];
            y = [crn(2,1) crn(2,1) crn(2,2) crn(2,2)];
            v = [x;y]';
            f = [1 2 3 4];
            if (x(1)<plotlim)
                rp = rp+1;
                x_rec(rp,1:4) = x;
                y_rec(rp,1:4) = y;
                %patch('Faces',f,'Vertices',v,'FaceColor','blue','FaceAlpha',0.2,'EdgeColor','None')
            end
            %text(x(1),(y(1)+y(3))/2,"E"+string(el),'Color','blue','FontSize',10)    
        end
        patch('XData',x_rec', 'YData',y_rec','FaceColor','blue','FaceAlpha',0.2,'EdgeColor','None')
    end
    
    %% create axis labels

    for i=1:length(mesh.ind_to_par{1})
        xt(i) = mesh.ind_to_par{1}(1,i); %#ok<AGROW>
        %xl(i) = string(mesh.ind_to_par{1}(1,i))+" \\ "+ string(mesh.ind_to_par{1}(2,i)); %#ok<AGROW>
        %xl(i) = "\begin{tabular}{c}"+ xl(i)+ "\end{tabular}"; %#ok<AGROW>
        if mesh.ind_to_par{1}(1,i)<9.5
            xl(i) =  string(mesh.ind_to_par{1}(1,i))+"\space\phantom{2}"+string(round(mesh.ind_to_par{1}(2,i),3)); %#ok<AGROW>
        else
            xl(i) =  string(mesh.ind_to_par{1}(1,i))+"\space"+string(round(mesh.ind_to_par{1}(2,i),3)); %#ok<AGROW>
        end
    end
    set(gca, 'TickLabelInterpreter', 'latex');
    ax = gca;
    xticks(xt)
    xticklabels(xl)
    xtickangle(-90)
    
    for i=1:length(mesh.ind_to_par{2})
        yt(i) = mesh.ind_to_par{2}(1,i); %#ok<AGROW>
        yl(i) = string(mesh.ind_to_par{2}(1,i)); %#ok<AGROW>
        for z=1:3-strlength(yl(i))
            yl(i) = yl(i)+"\phantom{2}";
        end
        st = string(round(mesh.ind_to_par{2}(2,i),3));
        for z=1:6-strlength(st)
            yl(i) = yl(i)+"\phantom{2}";
        end
        yl(i) = yl(i)+ string(round(mesh.ind_to_par{2}(2,i),3)); %#ok<AGROW>
        
    end
    yticks(yt)
    yticklabels(yl)
    
    ax.TickLength = [0.0 0.0];
    
    title("p="+string(mesh.p)+", parametric domain")
    
end

function plot_meshlines_even(mesh, plotlim)
    rpc = 0;

    for an=1:length(mesh.anchors)  
        knot_x = mesh.anchors{an}.knot_ind_x;
        knot_y = mesh.anchors{an}.knot_ind_y;
        
        i = floor(length(knot_x)/2);
        j = floor(length(knot_y)/2);

        corners = [knot_x(i), knot_y(j)];
        len = [knot_x(i+1)-knot_x(i), knot_y(j+1)-knot_y(j)];
        pos = [corners, len];
        
        if pos(1,1)<plotlim
            rpc = rpc+1;
            %rectangle('Position',pos,'EdgeColor','k')
            
            x_rec(rpc,1:4) = [knot_x(i) knot_x(i+1) knot_x(i+1), knot_x(i)];
            y_rec(rpc,1:4) = [knot_y(j) knot_y(j) knot_y(j+1) knot_y(j+1)];
        end
    end
    patch('XData',x_rec', 'YData',y_rec', 'EdgeColor','k','FaceColor','none')
    
end

function plot_meshlines_odd(mesh, plotlim)
    l=0;
    for an=1:length(mesh.anchors)
        knot_x = mesh.anchors{an}.knot_ind_x;

        
        %% left and right meshline
        if (knot_x(1)<plotlim+5)
            knot_y = mesh.anchors{an}.knot_ind_y;
        
            my_x = median(knot_x);
            my_y = median(knot_y);
            for cp=1:length(mesh.anchors)
                if (cp~=an)
                    cp_y = median(mesh.anchors{cp}.knot_ind_y);
                    if (cp_y == my_y) %points on same height
                        cp_knot_x = mesh.anchors{cp}.knot_ind_x;
                        if (ismember(my_x, cp_knot_x))
                            x = [my_x median(cp_knot_x)];
                            y = [my_y cp_y];
                            if (x(1)<plotlim)
                                l=l+1;
                                lx(l,1:2) = x;
                                ly(l,1:2) = y;
                                %line(x,y,'Color', 'k');
                            end
                        end
                    end
                end
            end

            %% up and down meshline
            for cp=1:length(mesh.anchors)
                if (cp~=an)
                    cp_x = median(mesh.anchors{cp}.knot_ind_x);
                    if (cp_x == my_x) %points on same height
                        cp_knot_y = mesh.anchors{cp}.knot_ind_y;
                        if (ismember(my_y, cp_knot_y))
                            x = [my_x cp_x];
                            y = [my_y median(cp_knot_y)];
                            if (x(1)<plotlim)
                                l=l+1;
                                lx(l,1:2) = x;
                                ly(l,1:2) = y;
                                %line(x,y,'Color', 'k');
                            end
                        end
                    end
                end
            end  
        end
    end
    
    line(lx',ly','Color','k')
end

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

function plot_basis(mesh, xplot, yplot)
    vals = zeros(length(xplot),length(yplot));
        
    for an=1:length(mesh.anchors)
        vals_spline = plot_global_basisfunc(mesh.anchors{an}.global_knots, xplot, yplot, true);
        hold on
        %plot3(anchors{an}.loc(1), anchors{an}.loc(2),max(max(vals)),'kx')
        vals = vals+vals_spline;
    end
    surf(xplot,yplot,vals')
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

function plot_bezier_extr(mesh, xplot, yplot)
    p = mesh.p;

    vals = zeros(length(xplot),length(yplot));

    for el=1:length(mesh.elementgroups{1}.elements)
        domain = mesh.elementgroups{1}.elements{el}.corners_par;
        Bez = mesh.elementgroups{1}.elements{el}.Bezier;
        
        B = get_bernstein(domain, p, xplot, yplot);
        for c=1:length(mesh.elementgroups{1}.elements{el}.cps)
            bez_extr = Bez(c,:);
            c_vals = zeros(length(xplot), length(yplot));
            for f=1:(p+1)^2
                c_vals = c_vals + B(:,:,f) * bez_extr(f);
            end
            [rows, cols] = find(c_vals);
            if (length(rows)<2 || length(cols)<2)
                
            else
                rows = [min(rows) max(rows)];
                cols = [min(cols) max(cols)];
                if (rows(1)~=rows(2) && cols(1)~=cols(2))
                    surf(xplot(rows(1):rows(2)),yplot(cols(1):cols(2)),c_vals(rows(1):rows(2), cols(1):cols(2))')
                end
            end
            hold on
            vals = vals+c_vals;
        end    
    end
    surf(xplot,yplot,vals')
    
end

function B = get_bernstein(domain, p, xplot, yplot)
    xknot = [domain(1,1) domain(1,2)];
    yknot = [domain(2,1) domain(2,2)];

    for i=1:p
        xknot = [xknot(1) xknot xknot(end)]; %#ok<AGROW>
        yknot = [yknot(1) yknot yknot(end)]; %#ok<AGROW>
    end
    
    Rx = zeros(p+1, length(xplot));
    Ry = zeros(p+1, length(yplot));
    for i=1:p+1
        knot_range(1) = i;
        knot_range(2) = knot_range(1)+p+1;
        
        b_knot_x = xknot(knot_range(1):knot_range(2));
        b_knot_y = yknot(knot_range(1):knot_range(2));
        
        b_spline_x = spmak(b_knot_x,1);
        b_spline_y = spmak(b_knot_y,1);

        Rx(i,:) = fnval(b_spline_x, xplot);
        Ry(i,:) = fnval(b_spline_y, yplot);
    end
    
    B = zeros(length(xplot),length(yplot),(p+1)^2);
    %checksum = zeros(length(xplot),length(yplot));
    for y=1:p+1
        for x=1:p+1
            B(:,:,x+(p+1)*(y-1)) = (Rx(x,:)'*Ry(y,:));
            %surf(xplot,yplot,(Rx(x,:)'*Ry(y,:))')
            hold on
            %checksum = checksum + (Rx(x,:)'*Ry(y,:))';
        end
    end
    %surf(xplot,yplot,checksum)
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
		ELS = [];
		DATA = [];
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

	h5disp(fname+".h5")
end
