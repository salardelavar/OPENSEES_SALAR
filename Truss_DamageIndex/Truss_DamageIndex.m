%%            IN THE NAME OF ALLAH           %%
%%        Truss Damage Index Analysis        %%
%%        Written by Salar Delavar Qashqai   %%
clc;close all;clear all;
t01=xlsread('Truss_DamageIndex_OUTPUT.csv','A2:A1001');%% Incerement
D03=xlsread('Truss_DamageIndex_OUTPUT.csv','F2:F1001');%% Disp. Incerement
DI01=xlsread('Truss_DamageIndex_OUTPUT.csv','AB2:AB1001');%% Damage Index element 01
DI02=xlsread('Truss_DamageIndex_OUTPUT.csv','AC2:AC1001');%% Damage Index element 02
DI03=xlsread('Truss_DamageIndex_OUTPUT.csv','AD2:AD1001');%% Damage Index element 03
DI04=xlsread('Truss_DamageIndex_OUTPUT.csv','AE2:AE1001');%% Damage Index element 04
DI05=xlsread('Truss_DamageIndex_OUTPUT.csv','AF2:AF1001');%% Damage Index element 05
DI06=xlsread('Truss_DamageIndex_OUTPUT.csv','AG2:AG1001');%% Damage Index element 02

Strain01=xlsread('Truss_DamageIndex_OUTPUT.csv','J2:J1001');%% Strain 01
Stress01=xlsread('Truss_DamageIndex_OUTPUT.csv','K2:K1001');%% Stress 01
Strain02=xlsread('Truss_DamageIndex_OUTPUT.csv','L2:L1001');%% Strain 02
Stress02=xlsread('Truss_DamageIndex_OUTPUT.csv','M2:M1001');%% Stress 02
Strain03=xlsread('Truss_DamageIndex_OUTPUT.csv','N2:N1001');%% Strain 03
Stress03=xlsread('Truss_DamageIndex_OUTPUT.csv','O2:O1001');%% Stress 03
Strain04=xlsread('Truss_DamageIndex_OUTPUT.csv','P2:P1001');%% Strain 04
Stress04=xlsread('Truss_DamageIndex_OUTPUT.csv','Q2:Q1001');%% Stress 04
Strain05=xlsread('Truss_DamageIndex_OUTPUT.csv','R2:R1001');%% Strain 05
Stress05=xlsread('Truss_DamageIndex_OUTPUT.csv','S2:S1001');%% Stress 05
Strain06=xlsread('Truss_DamageIndex_OUTPUT.csv','T2:T1001');%% Strain 06
Stress06=xlsread('Truss_DamageIndex_OUTPUT.csv','U2:U1001');%% Stress 06
figure(1)
plot(t01,DI01,t01,DI02,t01,DI03,t01,DI04,t01,DI05,t01,DI06,'LineWidth',2);
title('Damage Index Truss With 6 Elements During Pushover Analysis');
legend('ele.01','ele.02','ele.03','ele.04','ele.05','ele.06','Location','NorthEastOutside');
xlabel('Increment');ylabel('Damage Index');
figure(2)
%semilogy(D03,DI01,D03,DI02,D03,DI03,D03,DI04,D03,DI05,D03,DI06);
plot(D03,DI01,D03,DI02,D03,DI03,D03,DI04,D03,DI05,D03,DI06,'LineWidth',2);
title('Damage Index Truss With 6 Elements During Pushover Analysis');
legend('ele.01','ele.02','ele.03','ele.04','ele.05','ele.06','Location','NorthEastOutside');
xlabel('Displacement [DOF(3)]');ylabel('Damage Index');
figure(3)
plot(Strain01,Stress01,Strain02,Stress02,Strain03,Stress03,Strain04,Stress04,Strain05,Stress05,Strain06,Stress06,'LineWidth',2);
title('Strain-Stress of Truss With 6 Elements During Pushover Analysis');
legend('ele.01','ele.02','ele.03','ele.04','ele.05','ele.06','Location','NorthEastOutside');
xlabel('Strain');ylabel('Stress');