function [cal, cal_scaled] = calibrate_phase_position(exp, run, date, PV_name, PV_range)

    DAQ = [exp '_' run];
    data_path = ['/nas/nas-li20-pm00/' exp '/2024/2024' date '/' DAQ '/' DAQ '.mat'];
    
    load(data_path);
    
    R2 = data_struct.metadata.DTOTR2.RESOLUTION;
    comIndImg = data_struct.images.DTOTR2.common_index;
    comIndScal = data_struct.scalars.common_index;
    
    try
        PV_data = eval(['data_struct.scalars.' PV_name]);
    catch
        error('Specified PV does not exist in the data structure');
    end
    
    valid_indices = find(PV_data(comIndScal) >= PV_range(1) & PV_data(comIndScal) <= PV_range(2));
    
    xpos = [];
    phase = [];
    sigx = [];
    
    for i = valid_indices'
        img = double(imread(data_struct.images.DTOTR2.loc{comIndImg(i)}))';
        vec1 = sum(img, 1);
        [yfit, q, dq, chisq_ndf] = gauss_fit(1:length(vec1), vec1);
        
        try
            xpos(end + 1) = q(3) * R2;
            sigx(end + 1) = q(4) * R2;
            phase(end + 1) = data_struct.scalars.BSA_List_S20.TCAV_LI20_2400_P(comIndScal(i));
        catch
        end
    end
    
    f = 11.424e9; % Hz
    c = 3e8;      
    um_per_degRF = c / f / 360 * 1e6; 
    
    ifit = xpos < mean(xpos) * 2;
    p = polyfit(phase(ifit), xpos(ifit), 1);

    cal = abs(p(1));           % um/deg
    cal_scaled = cal / um_per_degRF; % um/um
    
end
