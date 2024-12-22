close all
clear
clc


commandXfoil = {'load k2.txt' ; 'ppar' ; 'n 364'; 'oper'};

ii = 5;
for a = -1.3:0.001:-1
    
    str1 = strcat('alfa', 32, num2str(a));
    str2 = 'cpwr';
    str3 = strcat('cp', num2str(a*1e4));
    commandXfoil(ii) = {str1};
    commandXfoil(ii+1) = {str2};
    commandXfoil(ii+2) = {str3};
    ii = ii + 3;
end

writecell(commandXfoil,'commandXfoil_cp.txt')









