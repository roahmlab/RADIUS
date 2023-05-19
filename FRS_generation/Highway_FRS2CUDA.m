clear, close all, clc

load('car_frs.mat');
m_mega = M_mega;
M_mega = cell(1,length(m_mega));


parfor i = 1:length(m_mega)
    M = containers.Map;
    M('Autb') = m_mega{i}('Autb');
    M('dirtb') = m_mega{i}('dirtb');
    M('lantb') = m_mega{i}('lantb');
    for mode = ["dir","Au","lan"]
        FRS = m_mega{i}(mode);
        tb = cell(size(FRS,1),size(FRS,2));
        if mode == "Au"
            manu_type = 1;
        elseif mode == "dir"
            manu_type = 2;
        else
            manu_type = 3;
        end

        for i1 = 1:size(FRS,1)
            for i2 = 1:size(FRS,2)
                
                frs = FRS{i1,i2};
                if isempty(frs)
                    continue
                end
                num_FRS_ini = length(frs.vehRS_save);

                % align time for the last FRS in vehRS_save
                Z = frs.vehRS_save{end}.Z;
                tc = Z(20,1);
                Z(20,2:end) = 0*Z(20,2:end);
                if abs(round(tc*100) - tc*100) > 1e-4
                    Z(20,2) = 0.005;
                else
                    Z(20,2) = 0.01;
                end
                frs.vehRS_save{end} = zonotope(Z);

                % start enclosing: try to enclose adjacent zonotope pair
                % into one zonotope
                frs_temp = {};
                cnt = 1;
                while cnt < length(frs.vehRS_save)
                    if cnt == length(frs.vehRS_save)
                        frs_temp{end+1} = frs.vehRS_save{end};
                    end
                    z = combine_zono(frs.vehRS_save{cnt}, frs.vehRS_save{cnt+1}, manu_type, 1);
                    intz = interval(z); intz = intz(1:2);
                    range_z = supremum(intz) - infimum(intz);
                    int1 = interval(frs.vehRS_save{cnt}); int1 = int1(1:2);
                    int2 = interval(frs.vehRS_save{cnt+1}); int2 = int2(1:2);
                    mini = min([infimum(int1),infimum(int2)],[],2);
                    maxi = max([supremum(int1),supremum(int2)],[],2);
                    range_ori = maxi - mini;
                    if all(range_z./range_ori <= 1.05)
                        frs_temp{end+1} = z;
                        cnt = cnt + 2;
                    else
                        frs_temp{end+1} = frs.vehRS_save{cnt};
                        cnt = cnt + 1;
                    end
                end
                frs.vehRS_save = frs_temp;

                % enclose one more time
                frs_temp = {};
                cnt = 1;
                while cnt < length(frs.vehRS_save)
                    if cnt == length(frs.vehRS_save)
                        frs_temp{end+1} = frs.vehRS_save{end};
                    end
                    z = combine_zono(frs.vehRS_save{cnt}, frs.vehRS_save{cnt+1}, manu_type);
                    intz = interval(z); intz = intz(1:2);
                    range_z = supremum(intz) - infimum(intz);
                    int1 = interval(frs.vehRS_save{cnt}); int1 = int1(1:2);
                    int2 = interval(frs.vehRS_save{cnt+1}); int2 = int2(1:2);
                    mini = min([infimum(int1),infimum(int2)],[],2);
                    maxi = max([supremum(int1),supremum(int2)],[],2);
                    range_ori = maxi - mini;
                    if all(range_z./range_ori <= 1.05)
                        frs_temp{end+1} = z;
                        cnt = cnt + 2;
                    else
                        frs_temp{end+1} = frs.vehRS_save{cnt};
                        cnt = cnt + 1;
                    end
                end
                frs.vehRS_save = frs_temp;

                % enclose one more time
                frs_temp = {};
                cnt = 1;
                while cnt < length(frs.vehRS_save)
                    if cnt == length(frs.vehRS_save)
                        frs_temp{end+1} = frs.vehRS_save{end};
                    end
                    z = combine_zono(frs.vehRS_save{cnt}, frs.vehRS_save{cnt+1}, manu_type);
                    intz = interval(z); intz = intz(1:2);
                    range_z = supremum(intz) - infimum(intz);
                    int1 = interval(frs.vehRS_save{cnt}); int1 = int1(1:2);
                    int2 = interval(frs.vehRS_save{cnt+1}); int2 = int2(1:2);
                    mini = min([infimum(int1),infimum(int2)],[],2);
                    maxi = max([supremum(int1),supremum(int2)],[],2);
                    range_ori = maxi - mini;
                    if all(range_z./range_ori <= 1.05)
                        frs_temp{end+1} = z;
                        cnt = cnt + 2;
                    else
                        frs_temp{end+1} = frs.vehRS_save{cnt};
                        cnt = cnt + 1;
                    end
                end
                frs.vehRS_save = frs_temp;



                frs.cuda_FRS = pre_process_FRS(frs.vehRS_save,manu_type,65);
                [i manu_type i1 i2 num_FRS_ini length(frs.vehRS_save)]

                tb{i1,i2} = frs;
            end
        end
        M(char(mode)) = tb;
    end
    M_mega{i} = M;
end

save('CUDA_Highway_frs.mat', 'M_mega','-v7.3','-nocompression')

