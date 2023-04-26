function singnalFiltered = LowPass(signal,fech, fLow,order)

[num,den]=butter(order,fLow/(fech/2),'low') ;
singnalFiltered = filtfilt(num,den,signal);