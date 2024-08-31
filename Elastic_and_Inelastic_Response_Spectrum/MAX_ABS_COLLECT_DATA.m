clc;
clear;
close all;
%% parameters
NumFiles=1000; % number of files
B = 200;% Width of Column Section
Ela = 200000;% Modulus of Elasticity of Section
Le = 3000;% Length of Column
A = B*B;% area of Section
I = (B*B*B*B)/12;% Secind Moment Intria of Section
Ke = (3*Ela*I)/power(Le,3);% Sitffness of Structure
NameFiles_D='DFree_'; % base name of files - DISPLACEMENT
NameFiles_V='VFree_'; % base name of files - VELOCITY
NameFiles_A='AFree_'; % base name of files - ACCELRATION
%% Code
MaxAbsCol_D=zeros(NumFiles,1);
MaxAbsCol_V=zeros(NumFiles,1);
MaxAbsCol_A=zeros(NumFiles,1);
PERIOD=zeros(NumFiles,1);
for i=1:NumFiles
Ma = 2*i;% Mass of Node in Structure    
PERIOD(i) = 2*3.1415*sqrt(Ma/Ke);% Sitffness of Structure    
filename_D = [NameFiles_D, mat2str(i), '.txt'];
A = textread(filename_D);
MaxAbsCol_D(i,1) = max(abs(A(:,2)));
filename_V = [NameFiles_V, mat2str(i), '.txt'];
A = textread(filename_V);
MaxAbsCol_V(i,1) = max(abs(A(:,2)));
filename_A = [NameFiles_A, mat2str(i), '.txt'];
A = textread(filename_A);
MaxAbsCol_A(i,1) = max(abs(A(:,2)));
end
%% Plot Data
figure(1)
p1=plot(PERIOD,MaxAbsCol_D);grid on;set(p1,'LineWidth',3);
%legend('PICK GROUND DISPLACEMENT','Nonlinear Analysis','Location','NorthEastOutside');
xlabel('PERIOD (t)');ylabel('PGD (mm)');
title(['MAX. PEACK GROUND DISPLACEMENT : ',num2str(max(MaxAbsCol_D)),' mm'],'color','b');
figure(2)
p1=plot(PERIOD,MaxAbsCol_V);grid on;set(p1,'LineWidth',3);
xlabel('PERIOD (t)');ylabel('PGV (mm)');
title(['MAX. PEACK GROUND VELOCITY : ',num2str(max(MaxAbsCol_V)),' mm/s'],'color','b');
figure(3)
p1=plot(PERIOD,MaxAbsCol_A);grid on;set(p1,'LineWidth',3);
xlabel('PERIOD (t)');ylabel('PGA (mm/s^2)');
title(['MAX. PEACK GROUND ACCELARATION : ',num2str(max(MaxAbsCol_A)),' mm/s^2'],'color','b');
%% Write Output
fileID=fopen('TOT_D.txt', 'w');
fprintf(fileID,'%f\n',MaxAbsCol_D);
fclose(fileID);
fileID=fopen('TOT_V.txt', 'w');
fprintf(fileID,'%f\n',MaxAbsCol_V);
fclose(fileID);
fileID=fopen('TOT_A.txt', 'w');
fprintf(fileID,'%f\n',MaxAbsCol_A);
fclose(fileID);


