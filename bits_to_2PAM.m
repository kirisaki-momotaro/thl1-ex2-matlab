function [ output ] = bits_to_2PAM( input )
output=input;

output( output == 1 )=-1;
output( output == 0 )=1;


end

