function Dates = findDates(InitialDate,step,time_interval,flag)
% Generate a date vector from an initial date given a time interval 
% in years and the step in months
%
% INPUT
% InitialDate:      starting date  
% step:             step in months     
% time_interval:    time interval in years
% flag:             0 -> include InitialDate in output vector 
%                   1 -> not include InitialDate in output vector
%
% OUTPUT
% Dates:            dates vector    
%
    switch flag
       case 0
           Dates = datetime(InitialDate,'ConvertFrom','datenum') ...
           + calmonths(0:step:time_interval*12);
           Dates = Dates';
       case 1
           Dates = datetime(InitialDate,'ConvertFrom','datenum') ...
           + calmonths(0:step:time_interval*12);
           Dates = Dates(2:end)';
    end
    
    % % Set the holidays
    % hols = holidays(Dates(1),Dates(end));
    % % Remove the President Day holiday
    % hols(hols == '21-Feb-2028') = [];
    % hols(hols == '21-Feb-2033') = [];
    % hols(hols == '20-Feb-2034') = [];
    % hols(hols == '21-Feb-2039') = [];
    % hols(hols == '20-Feb-2040') = [];
    % hols(hols == '20-Feb-2045') = [];
    % hols(hols == '21-Feb-2050') = [];
    % hols(hols == '20-Feb-2051') = [];
    % hols(hols == '21-Feb-2056') = [];
    % hols(hols == '21-Feb-2061') = [];
    % hols(hols == '20-Feb-2062') = [];
    % hols(hols == '21-Feb-2067') = [];
    % hols(hols == '20-Feb-2068') = [];

    % Business day check
    Dates=busdate(Dates-1,"follow",eurCalendar);

end %findDates