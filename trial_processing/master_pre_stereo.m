% processing code prior to stereo reconstructions
clear all
close all
clc

trial = 'TRM-09-06-2018-2334UTC';

% copy files after 
frames = 07200:07209;
copy_images(trial,frames)

% copy files after 
frames = 07200:11999;
copy_images(trial,frames)

% check if dropping images
check_image_drop(trial)