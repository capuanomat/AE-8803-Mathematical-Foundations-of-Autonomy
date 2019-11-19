%%%%% AE 8803: Autonomy - Homework 1 - Question 1 %%%%%
clear; clc;

% Boundary condition
t_f = 10;

% Initial condition
x_0 = [5, 2]';

% Dynamics and stuff
A = [0, 1; 0, 0];
B = [0; 1];
Q = [1, 0; 0, 1];
R = 1;
S = [2, 0; 0, 1];

% (a)
