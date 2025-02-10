module jlCCEtools

export get_decay_time

"""
    get_decay_time(time,signal)

Given a time/signal pair (both arrays must be equally long), return the 
time at which the initial signal has decayed to 1/e 
"""
function get_decay_time(time,signal)

    nsignal = size(signal)[1]
    ntime = size(time)[1]

    if nsignal != ntime
        error("get_decay_time: incompatible lengths for time and signal")
    end

    signal_div_e = signal[1]*exp(-1.)

    idx = -1
    for i in 1:nsignal
        if signal[i]<signal_div_e
            idx = i 
            break
        end
    end

    if idx == -1
        error("get_decay_time: cannot extrapolate, no decay to 1/e detected")
    end

    # get values at both sides of interval
    signal_0 = signal[idx-1]
    signal_1 = signal[idx]
    time_0 = time[idx-1]
    time_1 = time[idx]

    # linear coefficient
    m = (signal_1-signal_0)/(time_1-time_0)

    # get dt for improved timestep inside interval
    dt = (signal_div_e - signal_0) / m

    return time_0 + dt

end


end