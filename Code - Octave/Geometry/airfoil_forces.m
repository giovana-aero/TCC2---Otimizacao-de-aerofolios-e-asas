% Script utilizado pra fazer uma figura pro documento do TCC


clc,clear

%name = 'perfil_cst.txt';
name = 'naca2412.txt';
coo = dlmread(name);
%num = 41;

%ex = flip(coo(1:num,2));
%in = coo(num:end,2);
%camber = (z_c1-z_c2)/2;
%camber = (ex + in)/2;

figure(1),clf
n = 1;
s1=plot(coo(:,1),coo(:,2),'k','linewidth',n);,hold on,axis equal
%s2=plot([-0.2,1],[0,0],'b','linewidth',n);
s2=plot([-0,1],[0,0],'b','linewidth',n);
s5=plot([-0.2,0],[0,0]+0.002,'b--','linewidth',n);
s3=plot([-0.2,1],[-0.2,0],'k--','linewidth',n);
%s5=scatter(0.25,0,'k','filled');
%scatter((1-cosd(9.462322)),(1-sind(9.462322)),'k','filled')


circr = @(radius,rad_ang)  [radius*cos(rad_ang);  radius*sin(rad_ang)];         % Circle Function For Angles In Radians
circd = @(radius,deg_ang)  [radius*cosd(deg_ang);  radius*sind(deg_ang)];       % Circle Function For Angles In Degrees
N = 25;                                                         % Number Of Points In Complete Circle
r_angl = linspace(180, 180+9.462322, N);                             % Angle Defining Arc Segment (radians)
radius = 1.1;                                                   % Arc Radius
xy_r = circd(radius,r_angl); 
s4=plot(xy_r(1,:)+1, xy_r(2,:), 'k--') ;


direction = [0 0 1];
rotate(s1,direction,-9.462322)
rotate(s2,direction,-9.462322)
rotate(s3,direction,-9.462322)
rotate(s4,direction,-9.462322)
rotate(s5,direction,-9.462322)

%saveas(gcf,'revB_profile_forces_original.png')
x = cosd(9.462322)*0.75;
scatter((1-x),sind(9.65)*x-0.1,150,'k')



