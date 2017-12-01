%%test_run_linearized_mats.m
clear all
close all
clc
xin=[287,5,-176,0,2,0];
uin=[0,0];


[Al,Bl]=linearized_mats(xin,uin)