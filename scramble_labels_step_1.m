function [codebook , T_1] = scramble_labels_step_1(f , Pr , T_1 , numLevel , bit_index)
% Scramble the partition indexes. 
for u_index = 1 : length(T_1)
    primary_x = T_1(u_index , 2) ;
    binary_x = de2bi(primary_x - 1 , log2(numLevel) , 'left-msb') ;
    
    hold_var = binary_x(bit_index) ;
    binary_x(bit_index) = binary_x(1) ; 
    binary_x(1) = hold_var ; 
    
    secondary_x = bi2de(binary_x , 'left-msb') + 1;
    T_1(u_index , 2) = secondary_x ;
end
% Find the codebook based on the new partition indexes. Note that we only
% change/scramble the partition indexes, therefore, only the order of
% codewords in the codebook changes and the values of the codewords are not modified. 
parfor y = 1 : numLevel
    numerator = 0 ;
    denominator = 0 ;
    for x = 1 : numLevel
        u_index = find (T_1 (: , 2) == x) ;
        u = T_1(u_index , 1) ;
        
        numerator = numerator + Pr (x , y) * sum (u .* f(u_index)) ;
        denominator = denominator + Pr (x , y) * sum (f(u_index)) ;
    end
    codebook(y) = numerator / denominator ;
end
end