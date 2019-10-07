function compile_evaluate_elastic_energy

mex -v -O evaluate_elastic_energy.cpp -I/usr/include:/usr/local/include -L/usr/lib:/usr/local/lib -lCGAL -lgmp -lboost_thread

end
