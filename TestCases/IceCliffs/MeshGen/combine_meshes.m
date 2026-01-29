function [mesh, n_off1, n_off2] = combine_meshes(mesh1, mesh2, prefix1, prefix2)

   wb = waitbar(0,'Combining meshes');
    %% add prefix to elementgroup and nodegroup names for mesh1
    for i=1:length(mesh1.elementgroups)
        mesh.elementgroups{i} = mesh1.elementgroups{i};
        mesh.elementgroups{i}.name = prefix1+mesh.elementgroups{i}.name;
    end
    egroup_offset = i;
    
    for i=1:length(mesh1.nodegroups)
        mesh.nodegroups{i} = mesh1.nodegroups{i};
        mesh.nodegroups{i}.name = prefix1+mesh.nodegroups{i}.name;
    end
    ng_offset = i;

    %% no need to renumber mesh1
    mesh.anchors = mesh1.anchors;
    nodecount = length(mesh1.anchors);
    n_off1 = 0; n_off2 = nodecount;
    waitbar(0.1,wb,'Finished mesh 1')
    %% add prefix to elementgroup and nodegroup names for mesh2
    for i=1:length(mesh2.elementgroups)
        eg = i+egroup_offset;
        mesh.elementgroups{eg} = mesh2.elementgroups{i};
        mesh.elementgroups{eg}.name = prefix2+mesh.elementgroups{eg}.name;
    end
    for i=1:length(mesh2.nodegroups)
        ng = i+ng_offset;
        mesh.nodegroups{ng} = mesh2.nodegroups{i};
        mesh.nodegroups{ng}.name = prefix2+mesh.nodegroups{ng}.name;
    end
    waitbar(0.1,wb,'Added prefixes for mesh 2')
    %% renumber mesh2
	if (false)
	    for i=1:length(mesh2.anchors)
	        waitbar((i/length(mesh2.anchors)),wb,'Renumbering elements:  '+string(i)+'/'+string(length(mesh2.anchors)));
	        nodecount = nodecount+1;
	        an_old = i;
	        an_new = nodecount;
	        mesh.anchors{an_new} = mesh2.anchors{an_old};
	        for ng=1:length(mesh2.nodegroups)
	            ng_nodes = mesh2.nodegroups{ng}.nodes;
	            repl_locs = find(ng_nodes==an_old);
	            mesh.nodegroups{ng+ng_offset}.nodes(repl_locs) = an_new;     %#ok<FNDSB>
	        end
	        for eg=1:length(mesh2.elementgroups)
	            for el = 1:length(mesh2.elementgroups{eg}.elements)
	                eg_nodes = mesh2.elementgroups{eg}.elements{el}.cps;
	                repl_locs = find(eg_nodes==an_old);
	                mesh.elementgroups{eg+egroup_offset}.elements{el}.cps(repl_locs) = an_new;     %#ok<FNDSB>
	            end

	        end
	    end
	else 
    	old_to_new = zeros(2, length(mesh2.anchors));
    	for i=1:length(mesh2.anchors)
        	nodecount = nodecount+1;
        	an_old = i;
        	an_new = nodecount;
        	mesh.anchors{an_new} = mesh2.anchors{an_old};
        	old_to_new(1:2, i) = [an_old an_new];
    	end
    	waitbar(0.2,wb,'Renumbered cp for mesh 2')
	
    	for ng=1:length(mesh2.nodegroups)
        	ng_nodes = mesh2.nodegroups{ng}.nodes;
        	for i=1:length(ng_nodes)
            	ng_nodes(i) = old_to_new(2, old_to_new(1,:) == ng_nodes(i));
        	end
        	mesh.nodegroups{ng+ng_offset}.nodes = ng_nodes;
    	end
    	waitbar(0.3,wb,'Renumbered nodegroups for mesh 2')
    	
    	for eg=1:length(mesh2.elementgroups)
        	maxl = length(mesh2.elementgroups{eg}.elements);
        	for el = 1:maxl
            	eg_nodes = mesh2.elementgroups{eg}.elements{el}.cps;
            	for i=1:length(eg_nodes)
                	eg_nodes(i) = old_to_new(2, old_to_new(1,:) == eg_nodes(i));
            	end
            	mesh.elementgroups{eg+egroup_offset}.elements{el}.cps = eg_nodes;
            	if (mod(el, 100)==0)
                	waitbar(0.3+0.7*el/maxl ,wb,'Renumbering elements');
            	end
        	end
    	end
    
	end
    close(wb)
end