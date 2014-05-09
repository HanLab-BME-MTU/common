function output_strs = all_uncommon_str_takeout(input_strs)

stop_num = common_str_stop(input_strs);
cutback_num = common_str_cutback(input_strs);


for i_str = 1 : length(input_strs)
    output_strs{i_str} = input_strs{i_str}(stop_num:end-cutback_num);
end
