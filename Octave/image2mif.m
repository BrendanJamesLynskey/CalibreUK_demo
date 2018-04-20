% Script to write grayscale image as .mif file
% Prints unsigned hex
% Brendan Lynskey 2018

clear all;
close all;


im  = imread('lenna_256x256.bmp');
fid = fopen('lenna.mif','w'); 

imsize = size(im)
count = 0;

if (fid)     
    %% Write the RGB data
    fprintf(fid,'WIDTH = 8;\n');
    fprintf(fid,'DEPTH = %d;\n',imsize(1)*imsize(2));
    fprintf(fid,'ADDRESS_RADIX = HEX;\n');
    fprintf(fid,'DATA_RADIX = HEX;\n');
    fprintf(fid,'CONTENT BEGIN\n\n');
    for i=1:imsize(1)
       for j=1:imsize(2)
          fprintf(fid,'%x  : ',count);
          count = count + 1;
          fprintf(fid,'%x;\n', im(i,j));
       end    
    end 
    fprintf(fid,'END;\n');
    fclose(fid);
	
end



