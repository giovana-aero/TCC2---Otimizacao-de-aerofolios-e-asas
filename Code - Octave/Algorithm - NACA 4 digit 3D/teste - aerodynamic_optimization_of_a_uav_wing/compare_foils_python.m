clc,clear

foil1 = dlmread('foil_python.dat');
foil2 = dlmread('foil_octave.dat');
isequal(foil1,foil2)

%figure(1),clf
%plot(foil1(:,1),foil1(:,3),foil2(:,1),foil2(:,3)),grid on,axis equal