function [SDR_4 , T_4 , codebook ,Distortion] = ACOSQ_step_4(f , Pr , Pr_z, codebook , T_3 , numLevel , delta)
FileID = fopen ('Results.txt' , 'a') ;
D_4 = [2 1] ;
Threshold = 0.001 ;
while (D_4(1) - D_4(2)) / D_4(2) > Threshold/4
    D_4(1) = D_4(2) ;
    %% Optimal Partition 
    T_u = zeros(length(T_3) , 8) ;
    for u_index = 1 : length(T_3)
        summation = 0 ;
        d_4 = zeros(8 , 2) ;
        u = T_3(u_index , 1) ;
        for y_1 = 1 : 2
            for y_2 = 1 : 2
                for y_3 = 1 : 2
                    y_1_y_2 = (y_1 - 1) * 2 + y_2 ;
                    hold_x = T_3(u_index , 1 + y_1_y_2) ;
                    binary_x = de2bi(hold_x - 1 , log2(numLevel) , 'left-msb') ;
                    x_1 = binary_x(1) + 1 ;
                    x_2 = binary_x(2) + 1 ;
                    x_3 = binary_x(3) + 1 ;
                    x_1_x_2_x_3 = (x_1 - 1) * 4 + (x_2 - 1) * 2 + x_3 ;
                    
                    
                    y_1_y_2_y_3 = (y_1_y_2 - 1) * 2 + y_3 ;
                    for x_4 = 1 : 2
                        for y_4 = 1 : 2
                            y = (y_1_y_2_y_3 - 1) * 2 + y_4 ;
                            
                            summation = summation + ...
                                Pr_z(xor(x_3 - 1 , y_3 - 1) + 1 , xor(x_4 - 1 , y_4 - 1) + 1) *...
                                (u - codebook(y)) ^ 2 ;
                        end
                        d_4(y_1_y_2_y_3 , x_4) = summation ;
                        summation = 0 ;
                    end
                    [ ~ , partition_index] = min(d_4(y_1_y_2_y_3 , :)) ;
                    T_u (u_index , y_1_y_2_y_3 ) = (x_1_x_2_x_3 - 1) * 2 + partition_index ;
                end
            end
        end
    end
    T_4 = cat(2 , T_3(: , 1) , T_u) ;
    %% Optimal Centroids 
    for y_1_y_2_y_3 = 1 : 8
        for y_prime = 1 : 2
            numerator = 0 ;
            denominator = 0 ;
            y = (y_1_y_2_y_3 - 1) * 2 + y_prime ;
            for x = 1 : 16
                u_index = find(T_4(: , 1 + y_1_y_2_y_3) == x) ;
                u = T_4(u_index , 1) ;
                
                numerator = numerator + Pr(x , y) * sum(u .* f(u_index)) ;
                denominator = denominator + Pr(x , y) * sum(f(u_index)) ;
            end
            codebook(y) = numerator / denominator ;
        end
    end
    %% Distortion
    [D_4(2)] = distortion_4 (f , Pr , T_4 , codebook , numLevel , delta);
    fprintf (FileID , 'Overall D_4 = %f\n' ,D_4(2)) ;
end

SDR_4 = 10 * log10(1 / D_4(2)) ;
fprintf (FileID , 'Overall SDR_4 = %f\n' , SDR_4 ) ;
fprintf (FileID , '=================\n') ;
fclose (FileID) ;

Distortion = D_4(2) ;
end