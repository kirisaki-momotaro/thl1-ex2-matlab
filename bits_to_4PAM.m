function [ output ] = bits_to_4PAM( input )
output=input;

output( output == 1 )=-1;
output( output == 2 )=1;
output( output == 3 )=-3;
output( output == 4 )=3;


end

