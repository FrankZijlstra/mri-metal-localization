% Download all example datasets

clear all
close all

url = 'https://www.isi.uu.nl/People/Frank/mri-metal-localization/testdata.zip';
filename = './testdata.zip';

fprintf('Downloading test data...\n');
websave(filename, url);
fprintf('Unzipping test data...\n');
unzip(filename);
delete(filename);
fprintf('Done!\n');
