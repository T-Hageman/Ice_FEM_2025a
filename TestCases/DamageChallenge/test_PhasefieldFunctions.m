close all
clear all
clc

L = 2;
l = 0.1;
dx = l/10;
nx = round(L/dx)
meshOrder = 2;
DamageRegion = false;

Mesh = makeMesh(nx, L, meshOrder);
ncp = length(Mesh.spline.coefs);

for O=[4]
	if (O==1)
		clr = 'b';
	end
	if (O==2)
		clr = 'k';
	end
	if (O==4)
		clr = 'r';
	end

	P = zeros(ncp,1);
	for it=1:1
		[F,K,dist,int_dist] = getKF(Mesh, P, l, O, DamageRegion);
		
		du = - K\F;
		P = P + du;

		err = sum(abs(du.*F));
	end
	
	Mesh.spline.coefs = P;
	[F,K,dist,int_dist] = getKF(Mesh, P, l, O, DamageRegion);

	subplot(2,1,1)
	xplot = linspace(0,L,10*nx+1);
	yplot = fnval(Mesh.spline,xplot);
	plot(xplot, yplot, clr);
	hold on
	ylim([0 1.2])
	xlim([0 L])

	subplot(2,1,2)
	plot(Mesh.x_ip,dist, clr);
	hold on

	int_dist
end



function [F,K,dist, int_dist] = getKF(Mesh, P, l, O, DamageRegion)
	dummy = 10000.0;

	F = zeros(length(Mesh.spline.coefs),1);
	K = zeros(length(Mesh.spline.coefs));
	int_dist = 0;
	for i=1:length(Mesh.W)
		phi = Mesh.N(i,:)*P;
		gradphi = Mesh.G(i,:)*P;
		grad2phi = Mesh.G2(i,:)*P;

		if (O==1)
			if (phi>0)
				dist(i) = 3/8/l*(phi + l^2*gradphi^2);
				int_dist = int_dist + Mesh.W(i)*dist(i);
			else
				dist(i) = 0.0;
			end
	
			F = F + Mesh.W(i) * (  3/8/l* (Mesh.N(i,:)'*1 + 2*l^2*Mesh.G(i,:)'*gradphi)      );
			K = K + Mesh.W(i) * (  3/8/l* (                 2*l^2*Mesh.G(i,:)'*Mesh.G(i,:))  );
		end
		if (O==2)
			dist(i) = 0.5/l*(phi^2 + l^2*gradphi^2);
			int_dist = int_dist + Mesh.W(i)*dist(i);
	
			F = F + Mesh.W(i) * (  0.5/l*(2*Mesh.N(i,:)'*phi + 2*l^2*Mesh.G(i,:)'*gradphi)      );
			K = K + Mesh.W(i) * (  0.5/l*(2*Mesh.N(i,:)'*Mesh.N(i,:) + 2*l^2*Mesh.G(i,:)'*Mesh.G(i,:))  );
		end
		if (O==4)
			dist(i) = 1.0/19.2/l*(9.0*phi + 0.0*l^2*gradphi^2+l^4/1*grad2phi^2);
			int_dist = int_dist + Mesh.W(i)*dist(i);

			F = F + Mesh.W(i) * (  1.0/19.2/l*(Mesh.N(i,:)' + 0.0*l^2*Mesh.G(i,:)'*gradphi + 1/1*l^4*Mesh.G2(i,:)'*grad2phi)      );
			K = K + Mesh.W(i) * (  1.0/19.2/l*( 0.0*l^2*Mesh.G(i,:)'*Mesh.G(i,:) + 1/1*l^4*Mesh.G2(i,:)'*Mesh.G2(i,:))  );
		end

		if (abs(Mesh.x_ip(i)-Mesh.L/2)<Mesh.dx)
			F = F + dummy*Mesh.N(i,:)'*(1.0 - Mesh.N(i,:)*P);
			K = K - dummy*Mesh.N(i,:)'*Mesh.N(i,:);
		end
		if (abs(Mesh.x_ip(i)-0)<Mesh.dx)
			F = F + dummy*Mesh.N(i,:)'*(0.0 - Mesh.N(i,:)*P);
			K = K - dummy*Mesh.N(i,:)'*Mesh.N(i,:);
		end
		if (abs(Mesh.x_ip(i)-Mesh.L)<Mesh.dx)
			F = F + dummy*Mesh.N(i,:)'*(0.0 - Mesh.N(i,:)*P);
			K = K - dummy*Mesh.N(i,:)'*Mesh.N(i,:);
		end
	end

	%symmetry
	if (DamageRegion == false)
		F = F + dummy*Mesh.C_Centre'*(1.0 - Mesh.C_Centre*P);
		K = K - dummy*Mesh.C_Centre'*Mesh.C_Centre;
	end

end

function Mesh = makeMesh(nx, L, meshOrder)
	Mesh.L = L;
	Mesh.dx = L/nx;

	knts = [0,L];
	a = [1];
	for i=1:meshOrder
		knts = [0,knts,L];
		a = [a, 1];
	end

	nurbs = spmak(knts,a);
	ToAdd = linspace(0,L,nx+1); 
	ToAdd(ToAdd==L)=[];
	ToAdd(ToAdd==0)=[];
	nurbs = sprfn(nurbs,ToAdd);

	fnplt(nurbs);

	n_ip = meshOrder;
	[x1D, w1D] = getIpscheme(n_ip);

	xx_ip = [];

	for el=1:nx
		for i=1:n_ip
			x_ip = (el-1+x1D(i))*L/nx;
			xx_ip(end+1) = x_ip;
			w_ip = L/nx*w1D(i);

			for s=1:length(nurbs.coefs)
				nurbs.coefs = 0.0*nurbs.coefs;
				nurbs.coefs(s) = 1.0;

				N(n_ip*(el-1)+i,s) = fnval(nurbs, x_ip);
				G(n_ip*(el-1)+i,s) = fnval(fnder(nurbs,1), x_ip);
				G2(n_ip*(el-1)+i,s) = fnval(fnder(nurbs,2), x_ip);
			end
			W(n_ip*(el-1)+i) = w_ip;
		end
	end

	for s=1:length(nurbs.coefs)
		nurbs.coefs = 0.0*nurbs.coefs;
		nurbs.coefs(s) = 1.0;	
		C_Centre(s) = fnval(nurbs, 0.5*L);
	end

	Mesh.spline = nurbs;
	Mesh.x_ip = xx_ip;
	Mesh.N = N;
	Mesh.G = G;
	Mesh.G2=G2;
	Mesh.W = W;
	Mesh.C_Centre = C_Centre;
end

function [x1D, w1D] = getIpscheme(ipcount1D)
    if (ipcount1D == 1)
        x1D = 0;
        w1D = 2;
    elseif (ipcount1D == 2)
        x1D = [-1/sqrt(3); 1/sqrt(3)];
        w1D = [1; 1];                
    elseif (ipcount1D == 3)
        x1D = [-sqrt(3/5); 0; sqrt(3/5)];
        w1D = [5/9; 8/9; 5/9];                     
    elseif (ipcount1D ==4)
        x1D = [-sqrt(3/7+2/7*sqrt(6/5)); -sqrt(3/7-2/7*sqrt(6/5)); sqrt(3/7-2/7*sqrt(6/5)); sqrt(3/7+2/7*sqrt(6/5))];
        w1D = [(18-sqrt(30))/36; (18+sqrt(30))/36; (18+sqrt(30))/36; (18-sqrt(30))/36];                     
    elseif (ipcount1D == 5)
        x1D = [-1/3*sqrt(5+2*sqrt(10/7)); -1/3*sqrt(5-2*sqrt(10/7)); 0; 1/3*sqrt(5-2*sqrt(10/7)); 1/3*sqrt(5+2*sqrt(10/7))];
        w1D = [(322-13*sqrt(70))/900; (322+13*sqrt(70))/900; 128/225; (322+13*sqrt(70))/900; (322-13*sqrt(70))/900];                     
    else
        error("Higer order ip schemes not implemented in Shapes.Q9");
    end
    
    x1D = (x1D+1)/2;
	w1D = w1D/2;
end