function [x,y] = createProfile(Profilo,NPannelli,Chord)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here




        % Sfruttiamo XFoil per la creazione del profilo, sperando che sia
        % all'interno del database

        fileID = fopen('XFoilInput.txt','w');
        fprintf(fileID, ['naca ' ' '  Profilo, '\n\n']);
        fprintf(fileID,'pane\n\n');

        fprintf(fileID,'gdes\n');
        fprintf(fileID,'tgap 0 0 \n');
        
        fprintf(fileID,'exec \n\n\n');

        fprintf(fileID,'ppar\n');
        fprintf(fileID, ['n ' ' ' num2str(NPannelli+1)  '\n\n\n']);

        filename = strcat('NACA_', Profilo, '.dat');

        fprintf(fileID, ['save ' ' ' filename '\n\n']);
        fprintf(fileID,'y\n\n');
        fprintf(fileID,'quit \n\n');
        fclose(fileID);

        Str2Exec = strcat("xfoil < XFoilInput.txt > /dev/null 2>&1");
%         Str2Exec = strcat("xfoil < XFoilInput.txt ");

        system(Str2Exec);

        Corpo = importXfoilProfile(filename);

        % Prima flippa i vettori
        x = flipud(Corpo.x);
        y = flipud(Corpo.y);
        
        x = x.*Chord;
        y = y.*Chord;



end