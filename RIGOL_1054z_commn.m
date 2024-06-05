% myScope = oscilloscope;
% myScope.Resource = 'USB0::0x1AB1::0x04CE::DS1ZA252402124::INSTR';
% myScope.Driver = 'rgds1kz';
% connect(myScope);
% %%
% 
% %get(myScope) %To examine the current Oscilloscope setting
% 
% waveformArray = getWaveform(myScope, 'acquisition', true);
% waveformArray1 = readWaveform(myScope);
% 
% plot(waveformArray);% Plot the waveform.
% xlabel('Samples');
% ylabel('Voltage');
% 
% %%
% disconnect(myScope);
% clear myScope;
%% Code for communicating and analysing data from RIGOL DS1054 Z oscilloscope
clear();
clc;
% instrreset
oscScope = visa('ni', 'USB0::0x1AB1::0x04CE::DS1ZA252402124::INSTR');
oscScope.InputBufferSize = 1000000;
fopen(oscScope);

fprintf(oscScope, ':TRIG:SWE NORM'); % Set the oscilloscope to required trigger mode

maxVoltageValues = [];
maxVoltageValues2 = [];
maxVoltageFigure = figure(1);
ax = axes(maxVoltageFigure);
xlabel(ax, 'Shots');
ylabel(ax, 'Maximum Voltage (V)');
title(ax, 'Maximum Voltage vs. Shots');
hold(ax, 'on');
shots=0;
% Set up a loop to continuously acquire data after each trigger
while true
    % Check if there is a trigger event
    triggerStatus = query(oscScope, ':TRIGger:STATus?');
    
    if contains(triggerStatus, 'TD')
        % Read the waveform data from the oscilloscope
        fprintf(oscScope, ':WAV:SOUR CHAN1');  % Select the channel to read data from
        fprintf(oscScope, ':WAV:MODE NORM');   % Set the waveform data mode
        fprintf(oscScope, ':WAV:FORM ASC');
        % Read the waveform data and metadata
        waveformData = query(oscScope, ':WAV:DATA?');
        xIncrement = str2double(query(oscScope, ':WAV:XINC?'));
        xOrigin = str2double(query(oscScope, ':WAV:XOR?'));
        yIncrement = str2double(query(oscScope, ':WAV:YINC?'));
        yOrigin = str2double(query(oscScope, ':WAV:YOR?'));
%         
        fprintf(oscScope, ':WAV:SOUR CHAN2');  % Select the channel to read data from
        fprintf(oscScope, ':WAV:MODE NORM');   % Set the waveform data mode
        fprintf(oscScope, ':WAV:FORM ASC');
        % Read the waveform data and metadata
        waveformData2 = query(oscScope, ':WAV:DATA?');
        xIncrement2 = str2double(query(oscScope, ':WAV:XINC?'));
        xOrigin2 = str2double(query(oscScope, ':WAV:XOR?'));
        yIncrement2 = str2double(query(oscScope, ':WAV:YINC?'));
        yOrigin2 = str2double(query(oscScope, ':WAV:YOR?'));
        % Extract the numeric values from the waveform data
        %waveformValues = extractNumbers(waveformData(11:end)); % Skip the header
       %fprintf('Waveform data: %s\n', waveformData);
%         pattern = '-?\d+\.?\d*';
%         waveformValues = str2double(regexp(waveformData(11:end), pattern, 'match'));
        waveformValues = sscanf(waveformData(12:end),'%g,');
        waveformValues2 = sscanf(waveformData2(12:end),'%g,');
        % Calculate the time and voltage values
        timeValues = xOrigin + (0:length(waveformValues)-1) * xIncrement;
        voltageValues = yOrigin + waveformValues * yIncrement*12.5; %there is a scaling factor...Check again...
%         timeValues2 = xOrigin2 + (0:length(waveformValues2)-1) * xIncrement2;
        voltageValues2 = yOrigin2 + waveformValues2 * yIncrement2*12.5; %there is a scaling factor...Check again...
        % Find the maximum voltage value
        maxVoltage = max(voltageValues);
        maxVoltage2 = max(voltageValues2);
        % Append the maximum voltage value to the array
        maxVoltageValues = [maxVoltageValues, maxVoltage];
        maxVoltageValues2 = [maxVoltageValues2, maxVoltage2];
        shots=shots+1
        %plot maximum of data
        figure(1);
        plot(ax, shots, maxVoltage, 'bo', 'MarkerSize', 6);
        hold on
        plot(ax, shots, maxVoltage2, 'r*', 'MarkerSize', 6);
        drawnow;  
        % Plot the waveform data
        figure(2);
        plot(timeValues, voltageValues);
        hold on
        plot(timeValues, voltageValues2);
        title('Oscilloscope Waveform');
        xlabel('Time (s)');
        ylabel('Voltage (V)');
        drawnow;
    end

end
% Close the oscilloscope connection
fclose(oscScope);
delete(oscScope);
clear oscScope;
