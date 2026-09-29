# Here find the specific params used to run the pipelines

  + `realdata`: 
      + copy the params file for the analysis that want to run to `nf-realdata/params`
      (change output dir as wanted)
      + copy `../data/*/dataset/*` to `nf-realdata/input/`
  + `simulations`:
      + copy the params file for the simulations that want to run to `nf-simulations/params`
  + `simvc`:
      + copy the start_params that want (corresponding to the params file) to `nf-simulations/input/vc`
      + copy the h5 file in `../data/*/dataset/*.h5` to `nf-simulations/input/dataset/`
  
