clc
clear
close all

commandXfoil = {'load k2.txt' ; 'ppar' ; 'n 364';newline; 'oper'; 'visc' ; '100000'; 'iter 1000'};
ii=9;
r=[1e5 1e6];
a=[-1.21 0 2];
for p=1:length(r)
    str1 = strcat('r', 32, num2str(r(p)));    
    str2 = strcat('vpar');
    str3 = strcat('vacc 0.00001');
    commandXfoil(ii)= {str1};
    commandXfoil(ii+1)= {str2};
    commandXfoil(ii+2)= {str3};
    h=1;
    for n=5:4:13
        str4 = strcat('n',32,num2str(n));
        commandXfoil(ii+2+h) = {str2};
        commandXfoil(ii+2+h+1)= {str4};
        commandXfoil(ii+2+h+2)=  {' '};
        k=1;
        for f=1:length(a)
            str5 = strcat('alfa', 32, num2str(a(f)));
            %str6 = strcat('cpwr cpr', num2str(r(p)), 'a', num2str(a(f)*10),'n',num2str(n),'.txt');
            str7 = strcat('vplo');
            str8 = strcat('cf');
            str9 = strcat('dump cfr', num2str(r(p)), 'a', num2str(a(f)*100),'n',num2str(n),'.txt');
            str10 = strcat(' ');
            commandXfoil(ii+2+h+2+k) = {str5};
            %commandXfoil(ii+2+h+2+k+1) = {str6}; 
            commandXfoil(ii+2+h+2+k+1) = {str7};
            commandXfoil(ii+2+h+2+k+2) = {str8};
            commandXfoil(ii+2+h+2+k+3) = {str9};
            commandXfoil(ii+2+h+2+k+4) = {str10};
            k=k+5;
        end
        h=h+k+2;
    end
    ii=ii+h+2;
end
writecell(commandXfoil,'commandXfoil_grumman.txt')