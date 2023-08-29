function [new_path, evap_setting,num_path]=evap_setting_update(evap_setting,i,update_interval,trap_type,top_target,shot_mask)
% a function to adaptively update the evaporating sttings to keep a
% constant halo occupancy
if strcmp(trap_type,'quad_3_4_shunt_0_3')
    settings_list = {'c:\remote\settings202119Jul143339.xml';%evap: ~0.834 MHz
        'c:\remote\settings202119Jul141556.xml';%evap: ~0.836 MHz
        'c:\remote\settings202119Jul141352.xml';%evap: ~0.837 MHz
        'c:\remote\settings202119Jul141323.xml';%evap: ~0.838 MHz
        'c:\remote\settings202116Jul155525.xml';%evap: ~0.839 MHz
        'c:\remote\settings202113Jul095432.xml';%evap: ~0.840 MHz
        'c:\remote\settings202113Jul100347.xml';%evap: ~0.841 MHz
        'c:\remote\settings202119Jul151024.xml';%evap: ~0.841.5 MHz
        'c:\remote\settings202119Jul151147.xml';%evap: ~0.842 MHz
        'c:\remote\settings202119Jul151211.xml';%evap: ~0.842.5 MHz
        'c:\remote\settings202113Jul100553.xml';%evap: ~0.843 MHz
        'c:\remote\settings202119Jul153619.xml';%evap: ~0.843.5 MHz
        'c:\remote\settings202119Jul151256.xml';%evap: ~0.844 MHz
        'c:\remote\settings202113Jul101532.xml';%evap: ~0.845 MHz
        'c:\remote\settings202119Jul155149.xml';%evap: ~0.846 MHz
        };
elseif strcmp(trap_type,'quad_3_0_shunt_0_9')
%     settings_list = {
%         'c:\remote\settings202128Oct111018.xml';%evap: ~0.829 MHz
%         'c:\remote\settings202128Oct110923.xml';%evap: ~0.831 MHz
%         'c:\remote\settings202128Oct110843.xml';%evap: ~0.833 MHz
%         'c:\remote\settings202102Nov112644.xml';%evap: ~0.834 MHz
%         'c:\remote\settings202102Nov112450.xml';%evap: ~0.8345 MHz
%         'c:\remote\settings202128Oct110801.xml';%evap: ~0.835 MHz
%         'c:\remote\settings202102Nov112418.xml';%evap: ~0.8355 MHz
%         'c:\remote\settings202128Oct110721.xml';%evap: ~0.836 MHz
%         'c:\remote\settings202128Oct110636.xml';%evap: ~0.837 MHz
%         'c:\remote\settings202122Jul092216.xml';%evap: ~0.838 MHz
%         'c:\remote\settings202122Jul092456.xml';%evap: ~0.8385 MHz
%         'c:\remote\settings202120Jul151102.xml';%evap: ~0.839 MHz
%         'c:\remote\settings202122Jul092533.xml';%evap: ~0.8395 MHz
%         'c:\remote\settings202120Jul151039.xml';%evap: ~0.840 MHz
%         'c:\remote\settings202122Jul092648.xml';%evap: ~0.8405 MHz
%         'c:\remote\settings202120Jul151012.xml';%evap: ~0.841 MHz
%         'c:\remote\settings202120Jul150821.xml';%evap: ~0.8415 MHz
%         'c:\remote\settings202120Jul150601.xml';%evap: ~0.842 MHz
%         'c:\remote\settings202120Jul150519.xml';%evap: ~0.8425 MHz
%         'c:\remote\settings202120Jul150504.xml';%evap: ~0.843 MHz
%         'c:\remote\settings202120Jul150351.xml';%evap: ~0.8435 MHz
%         'c:\remote\settings202120Jul150327.xml';%evap: ~0.844 MHz
%         'c:\remote\settings202120Jul150310.xml';%evap: ~0.8445 MHz
%         'c:\remote\settings202120Jul142856.xml';%evap: ~0.845 MHz
%         'c:\remote\settings202120Jul151131.xml';%evap: ~0.846 MHz
%         'c:\remote\settings202120Jul151208.xml';%evap: ~0.847 MHz
%         'c:\remote\settings202122Jul161245.xml';%evap: ~0.848 MHz
%         'c:\remote\settings202122Jul161222.xml';%evap: ~0.849 MHz
%         'c:\remote\settings202122Jul161154.xml';%evap: ~0.850 MHz
%         'c:\remote\settings202122Jul161425.xml';%evap: ~0.851 MHz
%         'c:\remote\settings202122Jul161039.xml';%evap: ~0.852 MHz
%         'c:\remote\settings202122Jul161117.xml';%evap: ~0.853 MHz
%         'c:\remote\settings202122Jul161444.xml';%evap: ~0.854 MHz
%         };
    
% %updated in trapp cooling
%      settings_list = {
%         'c:\remote\settings202117Nov144103.xml';%evap: ~0.829 MHz
%         'c:\remote\settings202117Nov144131.xml';%evap: ~0.831 MHz
%         'c:\remote\settings202117Nov144150.xml';%evap: ~0.833 MHz
%         'c:\remote\settings202117Nov144208.xml';%evap: ~0.834 MHz
%         'c:\remote\settings202117Nov144257.xml';%evap: ~0.8345 MHz
%         'c:\remote\settings202117Nov144313.xml';%evap: ~0.835 MHz
%         'c:\remote\settings202117Nov144331.xml';%evap: ~0.8355 MHz
%         'c:\remote\settings202117Nov144355.xml';%evap: ~0.836 MHz
%         'c:\remote\settings202117Nov144416.xml';%evap: ~0.837 MHz
%         'c:\remote\settings202117Nov144439.xml';%evap: ~0.838 MHz
%         'c:\remote\settings202117Nov144529.xml';%evap: ~0.8385 MHz
%         'c:\remote\settings202117Nov144552.xml';%evap: ~0.839 MHz
%         'c:\remote\settings202117Nov144608.xml';%evap: ~0.8395 MHz
%         'c:\remote\settings202117Nov144627.xml';%evap: ~0.840 MHz
%         'c:\remote\settings202117Nov144649.xml';%evap: ~0.8405 MHz
%         'c:\remote\settings202117Nov144707.xml';%evap: ~0.841 MHz
%         'c:\remote\settings202117Nov144724.xml';%evap: ~0.8415 MHz
%         'c:\remote\settings202117Nov144801.xml';%evap: ~0.842 MHz
%         'c:\remote\settings202117Nov144817.xml';%evap: ~0.8425 MHz
%         'c:\remote\settings202117Nov144836.xml';%evap: ~0.843 MHz
%         'c:\remote\settings202117Nov144853.xml';%evap: ~0.8435 MHz
%         'c:\remote\settings202117Nov144910.xml';%evap: ~0.844 MHz
%         'c:\remote\settings202117Nov144925.xml';%evap: ~0.8445 MHz
%         'c:\remote\settings202117Nov144943.xml';%evap: ~0.845 MHz
%         'c:\remote\settings202117Nov145008.xml';%evap: ~0.846 MHz
%         'c:\remote\settings202117Nov145114.xml';%evap: ~0.847 MHz
%         'c:\remote\settings202117Nov145152.xml';%evap: ~0.848 MHz
%         'c:\remote\settings202117Nov145211.xml';%evap: ~0.849 MHz
%         'c:\remote\settings202117Nov145235.xml';%evap: ~0.850 MHz
%         'c:\remote\settings202117Nov145335.xml';%evap: ~0.851 MHz
%         'c:\remote\settings202117Nov145402.xml';%evap: ~0.852 MHz
%         'c:\remote\settings202117Nov145441.xml';%evap: ~0.853 MHz
%         'c:\remote\settings202117Nov145513.xml';%evap: ~0.854 MHz
%         };
%updated in trapp cooling (again)
%      settings_list = {
%         'c:\remote\settings202122Nov161747.xml';%evap: ~0.829 MHz
%         'c:\remote\settings202122Nov161734.xml';%evap: ~0.831 MHz
%         'c:\remote\settings202122Nov161701.xml';%evap: ~0.833 MHz
%         'c:\remote\settings202122Nov161254.xml';%evap: ~0.834 MHz
%         'c:\remote\settings202122Nov161814.xml';%evap: ~0.8345 MHz
%         'c:\remote\settings202122Nov161838.xml';%evap: ~0.835 MHz
%         'c:\remote\settings202122Nov161857.xml';%evap: ~0.8355 MHz
%         'c:\remote\settings202122Nov161916.xml';%evap: ~0.836 MHz
%         'c:\remote\settings202122Nov162004.xml';%evap: ~0.837 MHz
%         'c:\remote\settings202122Nov162019.xml';%evap: ~0.838 MHz
%         'c:\remote\settings202122Nov162142.xml';%evap: ~0.8385 MHz
%         'c:\remote\settings202122Nov162159.xml';%evap: ~0.839 MHz
%         'c:\remote\settings202122Nov162237.xml';%evap: ~0.8395 MHz
%         'c:\remote\settings202122Nov162256.xml';%evap: ~0.840 MHz
%         'c:\remote\settings202122Nov162314.xml';%evap: ~0.8405 MHz
%         'c:\remote\settings202122Nov162328.xml';%evap: ~0.841 MHz
%         'c:\remote\settings202122Nov162342.xml';%evap: ~0.8415 MHz
%         'c:\remote\settings202122Nov162357.xml';%evap: ~0.842 MHz
%         'c:\remote\settings202122Nov162414.xml';%evap: ~0.8425 MHz
%         'c:\remote\settings202122Nov162432.xml';%evap: ~0.843 MHz
%         'c:\remote\settings202122Nov162447.xml';%evap: ~0.8435 MHz
%         'c:\remote\settings202122Nov162506.xml';%evap: ~0.844 MHz
%         'c:\remote\settings202122Nov162536.xml';%evap: ~0.8445 MHz
%         'c:\remote\settings202122Nov162554.xml';%evap: ~0.845 MHz
%         'c:\remote\settings202122Nov162617.xml';%evap: ~0.846 MHz
%         'c:\remote\settings202122Nov162634.xml';%evap: ~0.847 MHz
%         'c:\remote\settings202122Nov162650.xml';%evap: ~0.848 MHz
%         'c:\remote\settings202122Nov162712.xml';%evap: ~0.849 MHz
%         'c:\remote\settings202122Nov162733.xml';%evap: ~0.850 MHz
%         'c:\remote\settings202122Nov162747.xml';%evap: ~0.851 MHz
%         'c:\remote\settings202122Nov162803.xml';%evap: ~0.852 MHz
%         'c:\remote\settings202122Nov162820.xml';%evap: ~0.853 MHz
%         'c:\remote\settings202122Nov162834.xml';%evap: ~0.854 MHz
%         };
     settings_list = {
        'c:\remote\settings202123Nov092834.xml';%evap: ~0.836 MHz
        'c:\remote\settings202123Nov092908.xml';%evap: ~0.837 MHz
        'c:\remote\settings202122Nov165903.xml';%evap: ~0.838 MHz
        'c:\remote\settings202122Nov165845.xml';%evap: ~0.8385 MHz
        'c:\remote\settings202122Nov165832.xml';%evap: ~0.839 MHz
        'c:\remote\settings202122Nov165819.xml';%evap: ~0.8395 MHz
        'c:\remote\settings202122Nov165801.xml';%evap: ~0.840 MHz
        'c:\remote\settings202122Nov165737.xml';%evap: ~0.8405 MHz
        'c:\remote\settings202122Nov165450.xml';%evap: ~0.841 MHz
        'c:\remote\settings202122Nov165717.xml';%evap: ~0.8415 MHz
        'c:\remote\settings202122Nov165605.xml';%evap: ~0.842 MHz
        'c:\remote\settings202122Nov165705.xml';%evap: ~0.8425 MHz
        'c:\remote\settings202122Nov165624.xml';%evap: ~0.843 MHz
        'c:\remote\settings202122Nov165645.xml';%evap: ~0.8435 MHz
        'c:\remote\settings202123Nov092753.xml';%evap: ~0.844 MHz
        };    
    
    
    
     settings_list_num = {
        'c:\remote\settings202102Aug161811.xml';%evap: ~0.838 MHz
        'c:\remote\settings202102Aug161757.xml';%evap: ~0.8385 MHz
        'c:\remote\settings202102Aug161743.xml';%evap: ~0.839 MHz
        'c:\remote\settings202102Aug161727.xml';%evap: ~0.8395 MHz
        'c:\remote\settings202102Aug161657.xml';%evap: ~0.840 MHz
        'c:\remote\settings202102Aug161642.xml';%evap: ~0.8405 MHz
        'c:\remote\settings202102Aug161609.xml';%evap: ~0.841 MHz
        'c:\remote\settings202102Aug161552.xml';%evap: ~0.8415 MHz
        'c:\remote\settings202102Aug161539.xml';%evap: ~0.842 MHz
        'c:\remote\settings202102Aug161523.xml';%evap: ~0.8425 MHz
        'c:\remote\settings202102Aug161449.xml';%evap: ~0.843 MHz
        'c:\remote\settings202102Aug161436.xml';%evap: ~0.8435 MHz
        'c:\remote\settings202102Aug161421.xml';%evap: ~0.844 MHz
        'c:\remote\settings202102Aug161405.xml';%evap: ~0.8445 MHz
        'c:\remote\settings202102Aug161323.xml';%evap: ~0.845 MHz
        'c:\remote\settings202102Aug161306.xml';%evap: ~0.846 MHz
        'c:\remote\settings202102Aug161251.xml';%evap: ~0.847 MHz
        'c:\remote\settings202102Aug161228.xml';%evap: ~0.848 MHz
        'c:\remote\settings202102Aug161154.xml';%evap: ~0.849 MHz
        'c:\remote\settings202102Aug161137.xml';%evap: ~0.850 MHz
        'c:\remote\settings202102Aug161113.xml';%evap: ~0.851 MHz
        'c:\remote\settings202102Aug161101.xml';%evap: ~0.852 MHz
        'c:\remote\settings202102Aug161030.xml';%evap: ~0.853 MHz
        'c:\remote\settings202102Aug160934.xml';%evap: ~0.854 MHz
        };
elseif strcmp(trap_type,'quad_0_7_shunt_0_75')
settings_list = {
    'c:\remote\settings202325Mar132858.xml';%832
    'c:\remote\settings202325Mar132842.xml';%834
    'c:\remote\settings202325Mar132817.xml';%836
    'c:\remote\settings202325Mar132747.xml';%838
    'c:\remote\settings202325Mar133039.xml';%839
    'c:\remote\settings202325Mar132731.xml';%840
    'c:\remote\settings202325Mar133001.xml';%841
    'c:\remote\settings202325Mar132711.xml';%evap: ~0.842 MHz
    'c:\remote\settings202325Mar132655.xml';%844
        'c:\remote\settings202316Mar164258.xml';%evap: ~0.846 MHz
        'c:\remote\settings202316Mar164154.xml';%evap: ~0.848 MHz
        'c:\remote\settings202316Mar164040.xml';%evap: ~0.850 MHz
        'c:\remote\settings202316Mar163923.xml';%evap: ~0.852 MHz
        'c:\remote\settings202316Mar163847.xml';%evap: ~0.854 MHz
        'c:\remote\settings202316Mar163814.xml';%evap: ~0.856 MHz
        'c:\remote\settings202316Mar163741.xml';%evap: ~0.858 MHz
        'c:\remote\settings202316Mar164339.xml';%evap: ~0.860 MHz
        'c:\remote\settings202318Mar214016.xml';%evap: ~0.861 MHz
        'c:\remote\settings202318Mar214130.xml';%evap: ~0.862 MHz
        'c:\remote\settings202318Mar214246.xml';%evap: ~0.863 MHz
        'c:\remote\settings202316Mar164404.xml';%evap: ~0.865 MHz
        'c:\remote\settings202316Mar164436.xml';%evap: ~0.870 MHz
        'c:\remote\settings202316Mar164459.xml';%evap: ~0.875 MHz
        'c:\remote\settings202316Mar170050.xml';%evap: ~0.880 MHz

        };
settings_list_num = {
        'c:\remote\settings202316Mar163814.xml';%evap: ~0.838 MHz
        };
end
if mod(i,update_interval) == update_interval-1
    %get the previous halos num
    
    if nargin>5
        halo_history = halo_num_check('Y:\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\',update_interval,shot_mask);
    else
        halo_history = halo_num_check('Y:\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\',update_interval);
    end
    num_check = (top_target - nanmean(halo_history.top.halo_N));
    if num_check<-1
        if num_check<-100
            evap_setting = evap_setting-10;
        elseif num_check<-18
            evap_setting = evap_setting-5;
        elseif num_check<-8
            evap_setting = evap_setting-3;
        else
            evap_setting = evap_setting-1;
        end
    elseif num_check>4
        if num_check>8
            evap_setting = evap_setting+2;
        else
            evap_setting = evap_setting+1;
        end
    end
end

if evap_setting<1
    evap_setting = 1;
elseif evap_setting>length(settings_list)
    evap_setting = length(settings_list);
end

new_path=settings_list{evap_setting};

if evap_setting>length(settings_list_num)
    num_path=settings_list_num{end};
else
    num_path=settings_list_num{evap_setting};
end


end