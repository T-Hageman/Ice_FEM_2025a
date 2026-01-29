function mesh = create_anchors_simple(p, repeats, NX, NY, surfaceMap)
    ny_per_layer = 4;
    ny_outer = 5;

    %% polynomial order
    mesh.p = p;
	mesh.Lx = 1;
	mesh.Ly = 1;
    Nx = NX;
    Ny = NY;
    
    %%  parameter mapping
    [ind_to_par_x, ind_to_par_y] = get_int_to_par(Nx, Ny, mesh.p, repeats);
    mesh.ind_to_par{1} = ind_to_par_x;
    mesh.ind_to_par{2} = ind_to_par_y;

    %%  knot vectors used for mesh
    [knotsx, knotsy] = get_knot_arrays(mesh, repeats, Nx, Ny);
    
    %%  locations for cp's
    [xlocs, ylocs] = cp_locations(mesh, knotsx, knotsy, repeats, 1, 1);

    %% create anchors
    an = 0;
    for j=1:length(knotsy{1})
        for i=1:length(knotsx{1})
            an=an+1;
            mesh.anchors{an}.knot_ind_x = knotsx{1}(i,:);
            mesh.anchors{an}.knot_ind_y = knotsy{1}(j,:); 
			xp = xlocs{1}(i)*(surfaceMap(1,end)-surfaceMap(1,1))+surfaceMap(1,1);
			yp = ylocs{1}(j)*interp1(surfaceMap(1,:),surfaceMap(2,:),xp);
            mesh.anchors{an}.cp = [xp, yp];
		end
	end
	mesh.Lx = max(surfaceMap(1,:));
	mesh.Ly = max(surfaceMap(2,:));
end

function [ind_to_par_x, ind_to_par_y] = get_int_to_par(Nx, Ny, p, repeats)
    %% get mapping from index to parametric
    %x direction
    loc = 0;
    ind_to_par_x = zeros(0,0);
    for i=1:floor(p/2)+1
        loc=loc+1;
        ind_to_par_x(:,loc) =[loc; 0];
    end
               
    for i=1:Nx-1                  
        for j=1:repeats
            loc=loc+1;
            ind_to_par_x(1,loc) = loc;
            ind_to_par_x(2,loc) = i/(Nx);
        end
    end
    
    for i=1:floor(p/2)+1
        loc=loc+1;
        ind_to_par_x(:,loc) =[loc; 1];
    end

    %y direction
    loc = 0;
    ind_to_par_y = zeros(0,0);
    for i=1:floor(p/2)+1
        loc=loc+1;
        ind_to_par_y(:,loc) =[loc; 0];
    end

    for i=1:Ny-1                  
        for j=1:repeats
            loc=loc+1;
            ind_to_par_y(1,loc) = loc;
            ind_to_par_y(2,loc) = i/(Ny);
        end
    end
    
    for i=1:floor(p/2)+1
        loc=loc+1;
        ind_to_par_y(:,loc) =[loc; 1];
    end

    
end

function [xlocs, ylocs] = cp_locations(mesh, knotsx, knotsy, repeats, LX, LY)
    p = mesh.p;
    
    %% x-direction
    for i=1:length(knotsx)
        uKnot = [0 0 1 1];
        i_add = unique(knotsx{i});
        i_add = i_add(mesh.p:length(i_add)-mesh.p+1);

        u_add = mesh.ind_to_par{1}(2,i_add);

		u_add = u_add(u_add~=0);
		u_add = u_add(u_add~=1);

        ControlP_in = zeros(2,2);
        ControlP_in(1,1) = 0.0;
        ControlP_in(1,2) = LX;

        Nb1 = nrbmak(ControlP_in,uKnot);
        Nb1 = nrbdegelev(Nb1,p-1);    
        Nb1 =  nrbkntins(Nb1,u_add);
        
        xlocs{i} = Nb1.coefs(1,:); %#ok<AGROW>
    end
    
    %% y-direction
    for i=1:length(knotsy)
        vKnot = [0 0 1 1];
        i_add = unique(knotsy{i});
        i_add = i_add(mesh.p:length(i_add)-mesh.p+1);

        v_add = mesh.ind_to_par{2}(2,i_add);
		v_add = v_add(v_add~=0);
		v_add = v_add(v_add~=1);


        ControlP_in = zeros(2,2);
        ControlP_in(1,1) = 0.0;
        ControlP_in(1,2) = LY;

        Nb2 = nrbmak(ControlP_in,vKnot);
        Nb2 = nrbdegelev(Nb2,p-1);    
        Nb2 =  nrbkntins(Nb2,v_add);
        
        ylocs{i} = Nb2.coefs(1,:); %#ok<AGROW>
    end

end

function  [knotsx, knotsy] = get_knot_arrays(mesh, repeats, Nx, Ny)

    ind_to_par_x = mesh.ind_to_par{1};
    ind_to_par_y = mesh.ind_to_par{2};

	fix = 0; fix2 = 0;
	if (mesh.p==3 && repeats == 2)
		%fix2 = 1;
		fix = 1;
	end
	if (mesh.p==2 && repeats == 1)
		fix = 1;
		%fix2 = 1;
	end
    
    p=mesh.p;
    %% x direction  
    n_elems_x = Nx;
    ymin = 1; 
    ymax = length(ind_to_par_x);
    for i=1:n_elems_x*repeats+p - (p-3) -fix
        for j=1:p+2
            knotsx{1}(i,j) = ind_to_par_x(1,min(max(j-ceil(p/2)-1+i-fix2,ymin),ymax));
        end
    end


    %% y direction  
    n_elems_y = Ny;
    ymin = 1; 
    ymax = length(ind_to_par_y);
    for i=1:n_elems_y*repeats+p - (p-3) -fix
        for j=1:p+2
            knotsy{1}(i,j) = ind_to_par_y(1,min(max(j-ceil(p/2)-1+i-fix2,ymin),ymax));
        end
    end

     
end