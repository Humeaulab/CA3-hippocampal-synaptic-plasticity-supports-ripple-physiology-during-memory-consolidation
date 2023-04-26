function singnalFiltered = BandPass(signal,fech, f1, f2,order)


singnalFiltered = HighPass(signal,fech, f1,order);
singnalFiltered = LowPass(singnalFiltered,fech, f2,order);