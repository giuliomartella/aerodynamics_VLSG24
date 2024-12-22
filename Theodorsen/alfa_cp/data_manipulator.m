close all
clear
clc

files = dir;

files = files(not([files.isdir]));
n=length(files(not([files.isdir])))-3;

A = importdata(files(1).name);
Points = str2double(A.textdata(4:end,1:2));
Angles = zeros(n,1);
Cp = zeros(length(A.textdata(4:end,1)),n);
for ii=1:n
    A = importdata(files(ii).name);
    Angles(ii) = str2double(A.textdata(2,3));
    Cp(:,ii) = str2double(A.textdata(4:end,3));
    ii
end


cpData.Points = Points;
cpData.Angles = Angles;
cpData.Cp = Cp;

save('cpData.mat',"cpData");


