clc,clear



%       iaf.designation = NACA 4 digit iaf.designation (eg. '2412') - STRING !
%                 iaf.n = no of panels (line elements) PER SIDE (upper/lower)
% iaf.HalfCosineSpacing = 1 for "half cosine x-spacing" 
%                       = 0 to give "uniform x-spacing"
%          iaf.wantFile = 1 for creating airfoil data file (eg. 'naca2412.dat')
%                       = 0 to suppress writing into a file
%       iaf.datFilePath = Path where the data  file has to be created
%                         (eg. 'af_data_folder/naca4digitAF/') 
%                         use only forward slash '/' (Just for OS portability)

% Data:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
%       af.x = x cordinate (nx1 array)
%       af.z = z cordinate (nx1 array)
%      af.xU = x cordinate of upper surface (nx1 array)
%      af.zU = z cordinate of upper surface (nx1 array)
%      af.xL = x cordinate of lower surface (nx1 array)
%      af.zL = z cordinate of lower surface (nx1 array)
%      af.xC = x cordinate of camber line (nx1 array)
%      af.zC = z cordinate of camber line (nx1 array)
%    af.name = Name of the airfoil
%  af.header = Airfoil name ; No of panels ; Type of spacing
%              (eg. 'NACA4412 : [50 panels,Uniform x-spacing]')

%iaf.designation = '2412';
%iaf.n = 40;
%iaf.HalfCosineSpacing = 1;
%iaf.wantFile = 0;
%iaf.datFilePath='./'; % Current folder
%iaf.is_finiteTE=0;
%af = naca4gen(iaf);
%
%figure(1),clf
%scatter(af.x,zeros(1,length(af.x))),hold on
%scatter(af.xU,zeros(1,length(af.xU))+0.1)
%
%coo = [];

x = cosspace_half(0,1,40);
m = 2;
p = 4; 
t = 12;
[xU,yU,xL,yL] = fourdigit(x,m,p,t);
coo = [flip(xU'),flip(yU');xL(2:end)',yL(2:end)'];
figure(1),clf
plot(coo(:,1),coo(:,2)),axis equal,grid on