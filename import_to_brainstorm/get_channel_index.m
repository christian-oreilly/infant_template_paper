function iChannel = get_channel_index(channels, channel_str)
    
    for iChannel=1:size(channels)
        
        if strcmp(channels(iChannel).Name, channel_str)
            return
        end
    end
    
    error(['Channel ' channel_str ' not found.'])
end
