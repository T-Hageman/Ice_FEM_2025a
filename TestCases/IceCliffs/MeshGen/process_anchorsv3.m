function [mesh] = process_anchorsv2(mesh)
   global wb;
   wb = waitbar(0,'Processing mesh');

   mesh = global_knots(mesh);
   mesh = create_nodegroups(mesh); 
   waitbar(0,wb,'Created global knots and nodegroups')

   mesh = create_interior_elements(mesh);
   mesh = Bezier_extract_interior(mesh);
   waitbar(0.99,wb,'Created interior elements') 
   
   mesh = create_boundary_elements(mesh);
   mesh = Bezier_extract_boundary(mesh);

   close(wb)
end

%% determine global knots
function mesh = global_knots(mesh) 
    %% global knot vectors in paramteric space
    lookup_x = mesh.ind_to_par{1};
    lookup_y = mesh.ind_to_par{2};
        
    for an=1:length(mesh.anchors)
        x_ind = mesh.anchors{an}.knot_ind_x;
        y_ind = mesh.anchors{an}.knot_ind_y;
        
        x_par = interp1(lookup_x(1,:), lookup_x(2,:), x_ind);
        y_par = interp1(lookup_y(1,:), lookup_y(2,:), y_ind);
        
        mesh.anchors{an}.global_knots = [x_par; y_par];
    end
end

%% create nodegroups
function mesh = create_nodegroups(mesh) 
    
    cp_loc = zeros(length(mesh.anchors),2);
    for an=1:length(mesh.anchors)
        cp_loc(an,1) = median(mesh.anchors{an}.knot_ind_x);
        cp_loc(an,2) = median(mesh.anchors{an}.knot_ind_y);  
    end
    
    min_x = min(cp_loc(:,1));
    min_y = min(cp_loc(:,2));
    max_x = max(cp_loc(:,1));
    max_y = max(cp_loc(:,2));
	mean_x = 0.5*(max_x-min_x);

    min_y_locs = find(cp_loc(:,2)==min_y);
    min_x_locs = find(cp_loc(:,1)==min_x);
    max_y_locs = find(cp_loc(:,2)==max_y);
    max_x_locs = find(cp_loc(:,1)==max_x);
    
    mesh.nodegroups{1}.name = "bottom";
    mesh.nodegroups{1}.nodes= sort(min_y_locs);
    mesh.nodegroups{2}.name = "right";
    mesh.nodegroups{2}.nodes= sort(max_x_locs);
    mesh.nodegroups{3}.name = "top";
    mesh.nodegroups{3}.nodes= sort(max_y_locs);
    mesh.nodegroups{4}.name = "left";
    mesh.nodegroups{4}.nodes= sort(min_x_locs);
    mesh.nodegroups{5}.name = "internal";
    mesh.nodegroups{5}.nodes= 1:1:length(mesh.anchors);
	mesh.nodegroups{6}.name = "sides";
    mesh.nodegroups{6}.nodes= [min_x_locs; max_x_locs];

	mesh.nodegroups{7}.name = "CentreBottom";
	mesh.nodegroups{7}.nodes = [];
	for i=1:length(cp_loc)
		if (abs(cp_loc(i,1)-mean_x)<=2 && cp_loc(i,2)==min_y)
			mesh.nodegroups{7}.nodes(end+1) = i;
		end
	end
end

%% create elements
function mesh = create_interior_elements(mesh)
    global wb;

    %% elements
    elgroup=1;
    mesh.elementgroups{elgroup}.name = "internal";
	mesh.elementgroups{elgroup}.type = "B2_" + string(length(mesh.anchors{1}.knot_ind_x)-2);
    tic;
    el = 0;
    elcorners_ind = zeros(0,0,0);
    for an=1:length(mesh.anchors)
        if (mod(an,25) == 0)
            s_toc = toc;
            progress = an/length(mesh.anchors);
            eta = s_toc *  (length(mesh.anchors) - an)/an;
            waitbar(progress,wb,'Creating interior elements:  '+string(an)+'/'+string(length(mesh.anchors)) + "   ETA:"+datestr(seconds(eta),'HH:MM:SS'));
        end
        
        xind  = mesh.anchors{an}.knot_ind_x;
        xknot = mesh.anchors{an}.global_knots(1,:);
        
        yind  = mesh.anchors{an}.knot_ind_y;
        yknot = mesh.anchors{an}.global_knots(2,:);
    
        for x=1:length(xind)-1
            for y=1:length(yind)-1
                add_element = false;
                new_element = false;
                duplicate_el = [];
                
                dx = xknot(x+1)-xknot(x);
                dy = yknot(y+1)-yknot(y);
                
                % check if element has a volume
                if (dx>0.0 && dy>0.0)
                    add_element = true;
                end
                
                % check first if element already exists
                if (add_element && el>0.0) 
                    duplicate = false;
                    duplicate_el = [];
					xind2 = [xind(x) xind(x+1)];
					yind2 = [yind(y) yind(y+1)];
                    for i=1:el %% parfor-red here
                        dx_ind = abs(xind2(1) - elcorners_ind(i,1,1) );
                        if (dx_ind<0.5)
                            dy_ind = abs(yind2(1) - elcorners_ind(i,2,1) );
                            if  (dy_ind<0.5)
                                duplicate_el = [duplicate_el, i]; 
                            end
                            dy_ind = abs(yind2(2) - elcorners_ind(i,2,2) );
                            if  (dy_ind<0.5)
                                duplicate_el = [duplicate_el, i]; 
                            end
                        end
                        
                        dx_ind = abs(xind2(2) - elcorners_ind(i,1,2) );
                        if (dx_ind<0.5)
                            dy_ind = abs(yind2(1) - elcorners_ind(i,2,1) );
                            if (dy_ind<0.5)
                                duplicate_el = [duplicate_el, i]; 
                            end
                            
                            dy_ind = abs(yind2(2) - elcorners_ind(i,2,2) );
                            if (dy_ind<0.5)
                                duplicate_el = [duplicate_el, i]; 
                            end
						end
						if (length(duplicate_el)==4)
							break;
						end
                    end
                    if (isempty(duplicate_el) == false)
                        duplicate = true;
                    end
                    
                    if (duplicate==false)
                        new_element = true;
                    end   
                elseif (add_element) %first element to be created, directly add-able
                    new_element = true;
                end
                
                % create new element
                if (new_element)
                    el=el+1;
                    mesh.elementgroups{elgroup}.elements{el}.corners_ind = [xind(x), xind(x+1);
                                                     yind(y), yind(y+1)];    
                    mesh.elementgroups{elgroup}.elements{el}.corners_par = [xknot(x), xknot(x+1);
                                                     yknot(y), yknot(y+1)];  
                    elcorners_ind(el,1:2,1:2) =  [xind(x), xind(x+1);
                                                  yind(y), yind(y+1)];                             
                    phys_1 = par_to_phys(mesh, [xknot(x), yknot(y)]);
                    phys_2 = par_to_phys(mesh, [xknot(x+1), yknot(y+1)]);
                                                 
                    mesh.elementgroups{elgroup}.elements{el}.corners_phys = [phys_1(1), phys_2(1);
                                                                             phys_1(2), phys_2(2)];
                    
                    mesh.elementgroups{elgroup}.elements{el}.cps = [];                         
                    duplicate_el = [el el el el];
                end
                
                % add anchor to relevant element
                if (add_element)
                    
                    if (length(unique(duplicate_el))==1 && length(duplicate_el)==4) %corresponds to a single (and complete) element
                       duplicate_el = unique(duplicate_el);
                       add_el = duplicate_el(1);
                       add_loc = length(mesh.elementgroups{elgroup}.elements{add_el}.cps)+1;
                       mesh.elementgroups{elgroup}.elements{add_el}.cps(add_loc) = an;
                    else % check if elements needs to be split
                       if (length(unique(duplicate_el))==1) %% element contained in large element, split element in parts if required
                           zz = duplicate_el(1);
                           dx_lb(1) = xind(x) - mesh.elementgroups{elgroup}.elements{zz}.corners_ind(1,1);
                           dx_lb(2) = yind(y) - mesh.elementgroups{elgroup}.elements{zz}.corners_ind(2,1);
                           
                           dx_rb(1) = xind(x+1) - mesh.elementgroups{elgroup}.elements{zz}.corners_ind(1,2);
                           dx_rb(2) = yind(y)   - mesh.elementgroups{elgroup}.elements{zz}.corners_ind(2,1);
                           
                           dx_lt(1) = xind(x)   - mesh.elementgroups{elgroup}.elements{zz}.corners_ind(1,1);
                           dx_lt(2) = yind(y+1) - mesh.elementgroups{elgroup}.elements{zz}.corners_ind(2,2);
                        
                           dx_rt(1) = xind(x+1) - mesh.elementgroups{elgroup}.elements{zz}.corners_ind(1,2);
                           dx_rt(2) = yind(y+1) - mesh.elementgroups{elgroup}.elements{zz}.corners_ind(2,2);
                           
                           if(dx_lb(1)<=0 && dx_rb(1)>=0 && dx_lb(2)<=0 && dx_lt(2)>=0) %pre-existing element is smaller
                               add_el = duplicate_el(1);
                               add_loc = length(mesh.elementgroups{elgroup}.elements{add_el}.cps)+1;
                               mesh.elementgroups{elgroup}.elements{add_el}.cps(add_loc) = an;
                           else %pre-exisitng element is larger
                               add_el = duplicate_el(1);
                               old_el = mesh.elementgroups{elgroup}.elements{add_el};
                               
                               add_loc = length(mesh.elementgroups{elgroup}.elements{add_el}.cps)+1;
                               mesh.elementgroups{elgroup}.elements{add_el}.cps(add_loc) = an;
                               
                               % duplicate element and re-determine corners
                               el=el+1;
                               mesh.elementgroups{elgroup}.elements{el} = old_el;
                               elcorners_ind(el,:,:) = elcorners_ind(add_el,:,:);
                               if (all(dx_lb == 0) && all(dx_lt == 0) ) %% split x direction
                                   mesh.elementgroups{elgroup}.elements{el    }.corners_ind(1,1) = xind(x+1);
                                   mesh.elementgroups{elgroup}.elements{add_el}.corners_ind(1,2) = xind(x+1);
                                   elcorners_ind(el,    1,1) = xind(x+1);
                                   elcorners_ind(add_el,1,2) = xind(x+1);
                                   mesh.elementgroups{elgroup}.elements{el    }.corners_par(1,1) = xknot(x+1);
                                   mesh.elementgroups{elgroup}.elements{add_el}.corners_par(1,2) = xknot(x+1);  
                                     
                                   phys = par_to_phys(mesh, [xknot(x+1), yknot(y+1)]);
                                                 
                                   mesh.elementgroups{elgroup}.elements{el    }.corners_phys(1,1) = phys(1);
                                   mesh.elementgroups{elgroup}.elements{add_el}.corners_phys(1,2) = phys(1);                                         
                               elseif (all(dx_rb == 0) && all(dx_rt == 0) )
                                   mesh.elementgroups{elgroup}.elements{el    }.corners_ind(1,2) = xind(x);
                                   mesh.elementgroups{elgroup}.elements{add_el}.corners_ind(1,1) = xind(x);
                                   elcorners_ind(el,    1,2) = xind(x);
                                   elcorners_ind(add_el,1,1) = xind(x);
                                   mesh.elementgroups{elgroup}.elements{el    }.corners_par(1,2) = xknot(x);
                                   mesh.elementgroups{elgroup}.elements{add_el}.corners_par(1,1) = xknot(x);  
                                
                                   phys = par_to_phys(mesh, [xknot(x), yknot(y)]);
                                                 
                                   mesh.elementgroups{elgroup}.elements{el    }.corners_phys(1,2) = phys(1);
                                   mesh.elementgroups{elgroup}.elements{add_el}.corners_phys(1,1) = phys(1);    
                                                               
                                   % check for zero span
                                   if (mesh.elementgroups{elgroup}.elements{el}.corners_par(1,2) == mesh.elementgroups{elgroup}.elements{el}.corners_par(1,1) )
                                       mesh.elementgroups{elgroup}.elements{el} = [];
                                       elcorners_ind(el,:,:) = [];
                                       el=el-1;
                                   end
                                   if (mesh.elementgroups{elgroup}.elements{add_el}.corners_par(1,2) == mesh.elementgroups{elgroup}.elements{add_el}.corners_par(1,1) )
                                       mesh.elementgroups{elgroup}.elements{add_el} = [];
                                       elcorners_ind(el,:,:) = [];
                                       el=el-1;
                                   end
                                   
                               else%% split y direction
                                   
							   end
                           end
                       elseif (length(duplicate_el))==4  %% all corners correspond to elements
                            unique_els = unique(duplicate_el);
                            for z=1:length(unique_els)
                               add_el = unique_els(z);
                               add_loc = length(mesh.elementgroups{elgroup}.elements{add_el}.cps)+1;
                               mesh.elementgroups{elgroup}.elements{add_el}.cps(add_loc) = an;             
							end
                       end  
                    end
                end
            end
        end
    end
    
    %reorder elements to be bottom-left top top-right numbered
    elements_old = mesh.elementgroups{elgroup}.elements;
    el = zeros(length(elements_old),1);
    for i=1:length(elements_old)
       el(i) = i;
       coords(i, :) =  elements_old{i}.corners_par(:,1);   
    end
    [~, re_ord] = sort(coords(:,2));
    el = el(re_ord);
    coords = coords(re_ord, :);
    
    rows = uniquetol(coords(:,2), 1e-6);
    ordered_els = [];
    for i=1:length(rows)
        locs = ismembertol(coords(:,2), rows(i), 1e-6);
        to_order_coords = coords(locs);
        to_order_el = el(locs);
        
        [~, re_ord] = sort(to_order_coords(:,1));
        ordered_els = [ordered_els; to_order_el(re_ord)]; %#ok<AGROW>
    end
    
    for i=1:length(ordered_els)
       el = ordered_els(i);
       mesh.elementgroups{elgroup}.elements{i} = elements_old{el}; 
    end
    
    %reorder to jive order (row-by-row, left-to-tight)
    for el=1:length(mesh.elementgroups{elgroup}.elements)
        cps = mesh.elementgroups{elgroup}.elements{el}.cps;
        an_reordered = cps*0;
        
        clear an coords
        for cp = 1:length(cps)
           an(cp) = cps(cp); %#ok<AGROW>
           coords(cp,:) = mesh.anchors{an(cp)}.cp;
        end
        [~, re_ord] = sort(coords(:,2));
        an = an(re_ord);
        coords = coords(re_ord, :);
        for i=1:mesh.p+1
            loc_ind = 1:1:mesh.p+1;
            loc_ind = loc_ind + (mesh.p+1)*(i-1);
            an_loc = an(loc_ind);
            coords_loc = coords(loc_ind,:);
            [~, re_ord] = sort(coords_loc(:,1));
            an_loc = an_loc(re_ord);
            for j=1:mesh.p+1
                an_reordered((i-1)*(mesh.p+1)+j) = an_loc(j);      
            end   
        end
        mesh.elementgroups{elgroup}.elements{el}.cps = an_reordered;
    end
    
end

function mesh = create_boundary_elements(mesh)
    
    for grp = 1:length(mesh.nodegroups)-3
        groupnum = length(mesh.elementgroups)+1;
        groupname = mesh.nodegroups{grp}.name;
        mesh.elementgroups{groupnum}.name = groupname;
		mesh.elementgroups{groupnum}.type = "B1_" + string(length(mesh.anchors{1}.knot_ind_x)-2);
        
        nodes = mesh.nodegroups{grp}.nodes;
        if (mesh.anchors{nodes(1)}.cp(2) == mesh.anchors{nodes(2)}.cp(2))
            horizontal = true;
        else
            horizontal = false;
        end
        mesh.elementgroups{groupnum}.horizontal = horizontal;
        
        el = 0;
        for nds = 1:length(nodes)
            an = nodes(nds);
            if (horizontal)
                ind   = mesh.anchors{an}.knot_ind_x;
                knots = mesh.anchors{an}.global_knots(1,:);
            else
                ind   = mesh.anchors{an}.knot_ind_y;
                knots = mesh.anchors{an}.global_knots(2,:);
            end
            for x=1:length(ind)-1
                add_element = false;
                new_element = false;
                
                dx = knots(x+1) - knots(x);
                
                if (dx>0) %% check if element has a length
                    add_element = true;
                end
                
                if (add_element == true)
                    if (el==0) %%directly create new element
                        new_element = true;
                    else  %% check if element is a duplicate
                        new_element = true;
                        for eln = 1:el
                            if (knots(x) == mesh.elementgroups{groupnum}.elements{eln}.corners_par(1) && ...
                                knots(x+1) == mesh.elementgroups{groupnum}.elements{eln}.corners_par(2))

                                duplicate = eln;
                                new_element = false;
                            end
                        end
                    end
                end
                
                if (new_element == true)
                    el=el+1;
                    
                    mesh.elementgroups{groupnum}.elements{el}.corners_ind  = [ind(x) ind(x+1)];
                    mesh.elementgroups{groupnum}.elements{el}.corners_par  = [knots(x) knots(x+1)];
                    %mesh.elementgroups{groupnum}.elements{el}.corners_phys = 
                    mesh.elementgroups{groupnum}.elements{el}.cps = [];
                 
                    duplicate = el;
                end
                
                if (add_element == true)
                   if (groupname == "top" || groupname == "left") %% account for direction rules of JIVE
                       mesh.elementgroups{groupnum}.elements{duplicate}.cps = [an mesh.elementgroups{groupnum}.elements{duplicate}.cps];
                   else
                       mesh.elementgroups{groupnum}.elements{duplicate}.cps(length(mesh.elementgroups{groupnum}.elements{duplicate}.cps)+1) = an;
                   end
                end 
            end       
        end       
	end   

	%% add combined side group
	groupnum = groupnum+1;
	mesh.elementgroups{groupnum}.name = "sides";
	mesh.elementgroups{groupnum}.type = "B1_" + string(length(mesh.anchors{1}.knot_ind_x)-2);
    mesh.elementgroups{groupnum}.horizontal = false;
	ec = 0;
	for g=[3,5]
		for el=1:length(mesh.elementgroups{g}.elements)
			ec = ec+1;
			mesh.elementgroups{groupnum}.elements{ec} = mesh.elementgroups{g}.elements{el};
		end
	end

end

function mesh = create_interface_elements(mesh)
    groupnum = length(mesh.elementgroups)+1;
    mesh.elementgroups{groupnum}.name = "interface";

    interface_index = find(mesh.ind_to_par{2}(2,:)==0.5);
    
    interface_nodes = [];
    for an=1:length(mesh.anchors)
        ind_y = mesh.anchors{an}.knot_ind_y;
        med = median(ind_y);
        if (med>interface_index(1+floor(mesh.p/2)-1) && med<interface_index(length(interface_index)-floor(mesh.p/2)+1))
            interface_nodes = [interface_nodes an]; %#ok<AGROW>
        end
    end

    el = 0;
    for nds = 1:length(interface_nodes)
            an = interface_nodes(nds);
            
            ind   = mesh.anchors{an}.knot_ind_x;
            knots = mesh.anchors{an}.global_knots(1,:);
            
            for x=1:length(ind)-1
                add_element = false;
                new_element = false;
                
                dx = knots(x+1) - knots(x);
                
                if (dx>0) %% check if element has a length
                    add_element = true;
                end
                
                if (add_element == true)
                    if (el==0) %%directly create new element
                        new_element = true;
                    else  %% check if element is a duplicate
                        new_element = true;
                        for eln = 1:el
                            if (knots(x) == mesh.elementgroups{groupnum}.elements{eln}.corners_par(1) && ...
                                knots(x+1) == mesh.elementgroups{groupnum}.elements{eln}.corners_par(2))

                                duplicate = eln;
                                new_element = false;
                            end
                        end
                    end
                end
                
                if (new_element == true)
                    el=el+1;
                    
                    mesh.elementgroups{groupnum}.elements{el}.corners_ind  = [ind(x) ind(x+1)];
                    mesh.elementgroups{groupnum}.elements{el}.corners_par  = [knots(x) knots(x+1)];
                    %mesh.elementgroups{groupnum}.elements{el}.corners_phys = 
                    mesh.elementgroups{groupnum}.elements{el}.cps = [];
                 
                    duplicate = el;
                end
                
                if (add_element == true)
                   mesh.elementgroups{groupnum}.elements{duplicate}.cps(length(mesh.elementgroups{groupnum}.elements{duplicate}.cps)+1) = an;
                end 
            end       
    end       

    for el = length(mesh.elementgroups{groupnum}.elements):-1:1
        if length(mesh.elementgroups{groupnum}.elements{el}.cps)== mesh.p+1 %non-discontinuous element
            mesh.elementgroups{groupnum}.elements(el) = [];                 % removing
        else %element contains discontinuity, append continuous nodes to obtain correct length
            newnodes_bot = mesh.elementgroups{groupnum}.elements{el}.cps(1:mesh.p+1);
            nodes_top = mesh.elementgroups{groupnum}.elements{el}.cps(mesh.p+2:end);
            for i=1:length(newnodes_bot)
                n_id = newnodes_bot(i);
                node_coords = mesh.anchors{n_id}.cp;
                n_top = -1;
                for j=1:length(nodes_top)
                   if (   node_coords(1) ==  mesh.anchors{nodes_top(j)}.cp(1)  )
                      n_top =  nodes_top(j);
                   end
                end
                if (n_top == -1)
                    newnodes_top(i) = n_id;
                else
                    newnodes_top(i) = n_top;
                end     
            end
            mesh.elementgroups{groupnum}.elements{el}.cps = [newnodes_bot newnodes_top];
        end
 
    end

end

%% Bezier extractor
function mesh = Bezier_extract_interior(mesh)
    %% bezier extractor for interior elements
    elgroup = 1;
    for el=1:length(mesh.elementgroups{elgroup}.elements)
        BE = zeros((mesh.p+1)^2, (mesh.p+1)^2); 
        
        cps = mesh.elementgroups{elgroup}.elements{el}.cps;  %% non-zero functions in element el
        par_corners = mesh.elementgroups{elgroup}.elements{el}.corners_par;
        interval_x = [par_corners(1,1) par_corners(1,2)];
        interval_y = [par_corners(2,1) par_corners(2,2)];
        %if (el==41)
        %   el
        %end
        for cp = 1:length(cps)
            anchor = cps(cp);
            knot = mesh.anchors{anchor}.global_knots;
            
            bez = bezier_extract(knot, interval_x, interval_y, mesh.p);
            
            BE(cp,1:(mesh.p+1)^2) = bez;
        end
        mesh.elementgroups{elgroup}.elements{el}.Bezier = BE;
    end

end

function mesh = Bezier_extract_boundary(mesh)

    for egroup = 2:length(mesh.elementgroups) % loop over boundary groups (1 is the interior)
        groupname = mesh.elementgroups{egroup}.name;
        for el = 1:length(mesh.elementgroups{egroup}.elements) % loop over all elements
            BE = zeros(mesh.p+1, mesh.p+1);
            cps = mesh.elementgroups{egroup}.elements{el}.cps;
            
            interval = mesh.elementgroups{egroup}.elements{el}.corners_par;
            for cp = 1:length(cps)
                an = cps(cp);
                if (mesh.elementgroups{egroup}.horizontal)
                    knots = mesh.anchors{an}.global_knots(1,:);
                else
                    knots = mesh.anchors{an}.global_knots(2,:);
                end
                %disp(groupname+': '+num2str(el)+'  '+num2str(an))
                bez = be_1D(knots, interval ,mesh.p);
                if (groupname == "top" || groupname == "left") %% account for direction rules of JIVE
                    BE(cp,1:mesh.p+1) = bez(end:-1:1);
                else
                    BE(cp,1:mesh.p+1) = bez;
                end
            end
            mesh.elementgroups{egroup}.elements{el}.Bezier = BE;
        end
        
        if (false) %for checking bezier extraction
            if (groupname == "top" || groupname == "left") %#ok<UNRCH>
                dir = -1;
            else
                dir = 1;
            end
            plot_bezier_1D(mesh.elementgroups{egroup}, mesh.p, dir); 
        end
    end

end

function mesh = Bezier_extract_interface(mesh)

    egroup = length(mesh.elementgroups); % interface is the last elementgroup
    for el = 1:length(mesh.elementgroups{egroup}.elements) % loop over all elements
        BE = zeros(mesh.p+1, mesh.p+1);
        cps = mesh.elementgroups{egroup}.elements{el}.cps;
            
        interval = mesh.elementgroups{egroup}.elements{el}.corners_par;
        for cp = 1:mesh.p+1
            an = cps(cp);
            knots = mesh.anchors{an}.global_knots(1,:);
            
            bez = be_1D(knots, interval ,mesh.p);
            BE(cp,1:mesh.p+1) = bez;
            
            mesh.elementgroups{egroup}.elements{el}.Bezier = BE;
        end
    end

end

function C = be_1D(knot, interval, p)
    

    if (interval(1)>=knot(1) && interval(2)<=knot(end) ) %interval is within knot vector
        el_index = find(unique(knot) == interval(2));
        el_index = el_index(1) - 1;

        knot_before = knot(knot~=knot(1));
        knot_before = knot_before(knot_before<interval(2));
        repeats = length(knot_before) - length(unique(knot_before));

        pre  = find(knot == knot(1));
        post = find(knot == knot(end));
        added = 0;

        while (length(pre)<p+1)
            added = added+1;
            knot = [knot(1) knot]; %#ok<AGROW>
            pre  = find(knot == knot(1));
        end
        while (length(post)<p+1)
            knot = [knot knot(end)]; %#ok<AGROW>
            post = find(knot == knot(end));
        end

        C = bezierExtraction(knot, p);

        func_numb = (added+1)-(el_index-1)-repeats;
        C = C(func_numb,:,el_index);
    else
        
       C = zeros(1,p+1);
        
        
    end

end

function bez = bezier_extract(knot, interval_x, interval_y, p)
    %% 1D bezier extractor on a per function, per interval basiS
    bez = zeros((p+1)^2,1);
	knot_x = knot(1,:);
	knot_y = knot(2,:);
    
    xmin = find(knot_x==interval_x(1));
    xmax = find(knot_x==interval_x(2));
    
    if (length(xmin)>=1 && length(xmax)>=1)
        C_x = be_1D(knot_x, interval_x, p);
        C_y = be_1D(knot_y, interval_y, p);

        
    else
        if (isempty(xmin))
           x_add = interval_x(1);
        else
           x_add = interval_x(2); 
        end
        
        Nb1 = spmak(knot_x, 1);
        Nb2 = fnrfn(Nb1, x_add);
        
        rat = Nb2.coefs;
        knot_app = Nb2.knots;
        
        C_x1 = be_1D(knot_app(1:p+2), interval_x, p);
        C_x2 = be_1D(knot_app(2:p+3), interval_x, p);
        C_x = C_x1.*rat(1) + C_x2.*rat(2);
        %spl = spmak(sort([interval_x interval_x interval_x interval_x]),C_x)
        %fnplt(spl)
        %hold on
        
        C_y = be_1D(knot_y, interval_y, p);
    end
    
    for j=1:p+1
        for i=1:p+1
           bez(i+(j-1)*(p+1)) = C_x(i)*C_y(j);
        end
    end
    
end

function plot_bezier_1D(egroup, p, dir) %#ok<DEFNU>
    xplot = linspace(0,1,762);
    vals = zeros(length(xplot),1);

    figure
    for el=1:length(egroup.elements)
        domain = egroup.elements{el}.corners_par;
        Bez = egroup.elements{el}.Bezier;
        
        B = get_bernstein_1D(domain, p, xplot);
        for c=1:length(egroup.elements{el}.cps)
            bez_extr = Bez(c,:);
            if (dir==-1)
               bez_extr = bez_extr(end:-1:1); 
            end
            c_vals = zeros(length(xplot), 1);
            for f=1:(p+1)
                c_vals = c_vals + B(:,f) * bez_extr(f);
            end
            plot(xplot,c_vals')
            hold on
            vals = vals+c_vals;
        end     
    end
    plot(xplot, vals')

end

%% related functions
function B = get_bernstein_1D(domain, p, xplot)
    xknot = [domain(1,1) domain(1,2)];

    for i=1:p
        xknot = [xknot(1) xknot xknot(end)]; %#ok<AGROW>
    end
    
    Rx = zeros(p+1, length(xplot));
    for i=1:p+1
        knot_range(1) = i;
        knot_range(2) = knot_range(1)+p+1; 
        b_knot_x = xknot(knot_range(1):knot_range(2));
        b_spline_x = spmak(b_knot_x,1);
        Rx(i,:) = fnval(b_spline_x, xplot);
    end
    
    B = zeros(length(xplot),(p+1));
    for x=1:p+1
        B(:,x) = Rx(x,:);
    end
end

function coords = par_to_phys(mesh, point)
    %% mapping from parametric to physical space
    
    if (true)
		%assuming perfectly rectangular
        coordsx = point(1)*mesh.Lx;
        coordsy = point(2)*mesh.Ly;
	else
        anchors = mesh.anchors;
        for an=1:length(mesh.anchors)
            global_knot = anchors{an}.global_knots;

            if (point(1)>=min(global_knot(1,:)) && point(1)<=max(global_knot(1,:)) && point(2)>=min(global_knot(2,:)) && point(2)<=max(global_knot(2,:)))
                cp = anchors{an}.cp;
                b_spline_x = spmak(global_knot(1,:),1);
                b_spline_y = spmak(global_knot(2,:),1);

                Rx = fnval(b_spline_x, point(1));
                Ry = fnval(b_spline_y, point(2));

                coordsx = coordsx+Rx*Ry*cp(1);
                coordsy = coordsy+Rx*Ry*cp(2);

                %"an:"+string(an)+" Rx: "+string(Rx)+" Ry: "+string(Ry)
            end
		end
    end
    coords = [coordsx coordsy];
end