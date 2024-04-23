
function [events_status,nsd_idx,shared_idx,nsd_repeats] = ieeg_nsdParseEvents(all_events)
%
% function [events_status,nsd_idx,shared_idx,nsd_repeats] = ieeg_nsdParseEvents(all_events)
% 
% Input
%       all_events: BIDS events table
%
% Output 
%   events_status 2 columns, 
%       column 1: 1 if the event is bad 
%       column 2: 1 for low interictal, 2 for high interictal
%
%   nsd_idx: 1530x1 column vector containing the nsd index numbers found
%   from the original file names for all shared images 
%   ('shared0814_nsd59080_prepped1.png' --> 59080). NaN: probe images.

%   shared_idx: 1530x1 column vector containing the shared index numbers
%   found from the original file names for all shared images
%   ('shared0814_nsd59080_prepped1.png' --> 814).  NaN: probe images.

%   nsd repeats: 0: nonrepeats, 1-6: repeats, NaN: probe
%
%
% DH, MM, Multimodal Neuroimaging Lab 2023


% events status (0:good, 1:bad) and interictal (0:no, 1:low, 2:high)
events_status = zeros(height(all_events),2);
    
% events status (0:good, 1:bad)
events_status(ismember(all_events.status,'bad'),1) = 1;

% find interictal label (0:no, 1:low, 2:high) from events_table status_description
for ii_evts = 1:height(all_events)
    if ~isempty(all_events.status_description{ii_evts})
        this_set_of_descriptions = strip(split(all_events.status_description{ii_evts},','));
        if sum(ismember(this_set_of_descriptions,'/Interictal-findings/Attribute/Categorical/Categorical-level/Low'))>0 
            events_status(ii_evts,2) = 1;
        elseif sum(ismember(this_set_of_descriptions,'/Interictal-findings/Attribute/Categorical/Categorical-level/High'))>0
            events_status(ii_evts,2) = 2;
        end
    end
end

% get NSD idx
nsd_idx = NaN(height(all_events),1);
shared_idx = NaN(height(all_events),1);
stim_image_file = all_events.stim_file;
temp = split(stim_image_file, '\');
temp = temp(:,end); % last split is filename
for ii_evts = 1:height(all_events)
    if ~isequal(temp{ii_evts}(1),'z') % note a probe trial
        temp_2 = extractBetween(temp{ii_evts},'nsd','_prepped');
        nsd_idx(ii_evts) = str2double(temp_2{1});
        temp_3 = extractBetween(temp{ii_evts},'shared','_nsd');
        shared_idx(ii_evts) = str2double(temp_3{1});
    end 
end
clear temp_2 temp_3

% code repeats, 0: nonrepeats, 1-6: repeats, NaN: probe
nsd_repeats = NaN(height(all_events),1);
for ii_evts = 1:height(all_events)
    if length(find(nsd_idx == nsd_idx(ii_evts)))>1 % repeat
        nsd_repeats(ii_evts) = length(find(nsd_idx(1:ii_evts) == nsd_idx(ii_evts))); % how many times has it occured before ii_evts
    elseif length(find(nsd_idx == nsd_idx(ii_evts)))==1 % nonrepeat
        nsd_repeats(ii_evts) =  0;
    end
end
