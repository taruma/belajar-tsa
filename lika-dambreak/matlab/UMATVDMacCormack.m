%% INFORMATION! READ ME FIRST
% Program Saint-Venant: TVD-MacCormack Scheme
% Case: Dam Break
% Script by Taruma (25017046), hi@taruma.info
% Documentation: N/A
% Github this project: https://github.com/taruma/belajar-tsa
% Folder: lika-dambreak/matlab
% WARNING, this script is using Function in script, Use 2016b or later!!

%% --- Initiating Program
close all; clc; clear;

%% --- Define Constant
g = 9.81; beta = 1; n = 0.03;

%% --- Channel
h1 = 1.5; h2 = 1; lf = 55; m = 0; b = 5;

%% --- Secant
yawal = 5; dy = 0.1; yit = 10;

%% --- Grid and Time
dx = 1;
dt = 0.1;
time = 50;
tone = 1/dt; tout = tone*1;
imax = lf/dx+1;
tmax = time/dt;

