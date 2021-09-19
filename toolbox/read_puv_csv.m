function puv = read_puv_csv(puvfile)

PUV         = readtable(puvfile);
puv.freq    = PUV{:,1};
puv.SSE     = PUV{:,2};
puv.th      = PUV{:,3};
puv.sig     = PUV{:,4}*180/pi;