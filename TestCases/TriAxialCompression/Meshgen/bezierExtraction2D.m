function C = bezierExtraction2D(uknot,vknot,p,q)
	%
	% Bezier extraction operators for a 2D NURBS/BSpline.
	%
	% VP Nguyen
	% Cardiff University, UK
	
	% Bezier extraction operators for xi and eta
	% nb1: number of elements along xi direction
	% nb2: number of elements along eta direction
	
	[Cx,nb1]  = bezierExtraction(uknot,p);
	[Cy,nb2]  = bezierExtraction(vknot,q);
	
	% For Bsplines/NURBS, the element Bezier extraction
	% operator is square.

	% Bezier extraction operators for the whole mesh
	% as the tensor product of Cxi and Cet

	C = cell(nb1,nb2);
	for y=1:nb2
		for x=1:nb1
			BE = kron(Cy(:,:,y),Cx(:,:,x));
			if (~isempty(BE))
				C{x,y}(:,:) = BE;
			end
		end
	end
end