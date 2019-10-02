function [SDR_3 , T_3 , codebook] = ACOSQ_step_3(f , Pr , Pr_z , codebook , T_2 , numLevel , delta)
FileID = fopen ('Results.txt' , 'a') ;
Threshold = 0.001 ;
D_3 = [2 1] ;
while (D_3(1) - D_3(2)) / D_3(2) > Threshold/4
    D_3(1) = D_3(2) ;
    %% Optimal Partitions
    T_u = zeros(length(T_2) , 4) ;
    for u_index = 1 : length(T_2)
        d_3 = zeros(4 , 4) ;
        summation = 0 ;
        u = T_2(u_index , 1) ;
        for y_1 = 1 : 2
            
            hold_x = T_2(u_index , 1 + y_1) ;
            
            binary_x = de2bi(hold_x - 1 , log2(numLevel) , 'left-msb');
            x_1 = binary_x(1) + 1;
            x_2 = binary_x(2) + 1;
            x_1_x_2 = (x_1 - 1) * 2 + x_2 ;
            for y_2 = 1 : 2
                y_1_y_2 = (y_1 - 1) * 2 + y_2 ;
                for x_3 = 1 : 2
                    for x_4 = 1 : 2
                        x_prime = (x_3 - 1) * 2 + x_4 ;
                        for y_3 = 1 : 2
                            for y_4 = 1 : 2
                                
                                y = (y_1_y_2 - 1) * 4 + (y_3 - 1) * 2 + y_4 ;
                                summation = summation + ...
                                    Pr_z(xor(x_2 - 1 , y_2 - 1) + 1 , xor(x_3 - 1 , y_3 - 1) + 1) * ...
                                    Pr_z(xor(x_3 - 1 , y_3 - 1) + 1 , xor(x_4 - 1 , y_4 - 1) + 1) * ...
                                    (u - codebook(y)) ^ 2 ;
                            end
                        end
                        d_3(y_1_y_2 , x_prime) = summation ;
                        summation = 0 ;
                    end
                end
                [~ , partition_index] = min(d_3(y_1_y_2 , :)) ;
                T_u(u_index , y_1_y_2 ) = (x_1_x_2 - 1) * 4 + partition_index ;
            end
        end
    end
    
    T_3 = cat(2 , T_2(: , 1) , T_u) ;
    %% Optimal Centroids 
    for y_1_y_2 = 1 : 4
        for y_prime = 1 : 4
            y = (y_1_y_2 - 1) * 4 + y_prime ;
            numerator = 0 ;
            denominator = 0 ;
            for x = 1 : numLevel
                u_index = find(T_3(: , 1 + y_1_y_2) == x) ;
                u = T_3(u_index , 1) ;
                
                numerator = numerator + Pr(x , y) * sum(u .* f(u_index)) ;
                denominator = denominator + Pr(x , y) * sum(f(u_index)) ;
                
            end
            codebook(y) = numerator / denominator ;
        end
    end
    %% Distortion
    [D_3(2)] = distortion_3(f , Pr , T_3 , codebook , numLevel , delta) ;
    fprintf (FileID , 'Overall D_3 = %f\n' ,D_3(2)) ;
end

fprintf (FileID , 'Overall D_3 = %f\n' ,D_3(2)) ;
SDR_3 = 10 * log10(1 / D_3(2)) ;
fprintf (FileID , 'Overall SDR_3 = %f\n' , SDR_3 ) ;

fprintf (FileID , '=================\n') ;
fclose (FileID) ;
end