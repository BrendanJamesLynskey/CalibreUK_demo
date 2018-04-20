% Script to write signal array as .mif file
% Prints signed decimals
% Brendan Lynskey 2018

clear all;
close all;


sig_len   = 256;
sig_bits  = 8;
sig_max   = 2^sig_bits - 1;

sig_range = 0:(sig_len-1);
sig_range /= sig_len;
sig       = sig_max * sin(2.*pi.*sig_range);
sig       = round(sig);

fid       = fopen('sin.mif','w'); 


size(sig)

count     = 0;

if (fid)     
    %% Write the RGB data
    fprintf(fid,'WIDTH = %d;\n', sig_bits);
    fprintf(fid,'DEPTH = %d;\n', sig_len);
    fprintf(fid,'ADDRESS_RADIX = HEX;\n');
    fprintf(fid,'DATA_RADIX = DEC;\n');
    fprintf(fid,'CONTENT BEGIN\n\n');
    for i=1:size(sig)(2)
        fprintf(fid,'%x  : ',count);
        count = count + 1;
        fprintf(fid,'%d;\n', sig(1,i));
    end 
    fprintf(fid,'END;\n');
    fclose(fid);
	
end



