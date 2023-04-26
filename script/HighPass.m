function singnalFiltered = HighPass(signal,fech, fHigh,order)

[num,den]=butter(order,fHigh/(fech/2),'high') ;
singnalFiltered = filtfilt(num,den,signal);