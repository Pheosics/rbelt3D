	clear all;
	format long;

	box=6;

	xmin=-box;
	xmax=box;
	ymin=-box;
	ymax=box;
	zmin=-box;
	zmax=box;

	%basename='TS05-20031026';
	basename='dipole';
	fnumber=0;

	filename=[basename,'-dist-',num2str(fnumber),'.dat']

	fp = fopen(filename,'r');
	for j = 1:rbelt_info_read(basename,fnumber)
		[line,c] = fscanf(fp, '%u %u %e %e %e %e %e %e %e %e %e\n', [11,1]);
		x1(j) = line(3,1);
		y1(j) = line(4,1);
		z1(j) = line(5,1);
	end
	fclose(fp);

	f1=figure;
	plot3(x1,y1,z1,'k','LineWidth',0.2);
	%hold on;
	%plot3(x1,y1,z1,'*r');

%	hold on;
%	plot3(0.5000E+01,  0.1762E-01,  0.6838E-01,'*r');
%	hold on;
%	plot3(0.4999E+01,  0.1767E-01,  0.6873E-01,'*r');

	hold on;
	[xx,yy,zz]  = sphere(25);
	Cm=bone(50);
	surf(xx,yy,zz,xx,'EdgeColor','w')
	mesh(xx,yy,zz,xx)
	colormap(Cm);
	hold off;
	xlabel('x (R_E)','FontSize',12);
	ylabel('y (R_E)','FontSize',12);
	zlabel('z (R_E)','FontSize',12);
	axis equal;
	grid on;

	axis([xmin xmax ymin ymax zmin zmax]);
	view([1 -1 .4]);
	%view([0 0 1]);

	%axis([-5.5 -4.5 -1 1 -1 1]);
	%view([1 0 0]);



