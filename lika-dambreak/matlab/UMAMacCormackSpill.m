%% INFORMATION! READ ME FIRST
% Program Saint-Venant: MacCormack Scheme
% Case: .. - Spillway
% Script by Taruma (25017046), hi@taruma.info
% Documentation: N/A
% Github this project: https://github.com/taruma/belajar-tsa
% Folder: lika-dambreak/matlab

close all;
clear;

%% UH Q-time (m3/s-minute)

UH_Q = [0, 100, 0];
UH_time = [0, 500, 2500];
grad_Q = zeros(1, length(UH_time));

for i = 1:length(UH_Q)
    if i < length(UH_Q)
        grad_Q(i) = (UH_Q(i+1)-UH_Q(i))/(UH_time(i+1)/UH_time(i));
    else
        grad_Q(i) = 0;
    end
end