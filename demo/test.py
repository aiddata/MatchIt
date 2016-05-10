from mpi4py import MPI
from subprocess import Popen, PIPE, STDOUT, call, check_output
import random

comm = MPI.COMM_WORLD

size = comm.Get_size()
rank = comm.Get_rank()


command = 'Rscript'

sim_path = '/sciclone/home00/geogdan/MatchIt/demo/pySims.R'

iterations = 100
c = rank
while c < iterations:
    out_path = '/sciclone/home00/geogdan/AlphaSims/test_'+str(c)+'.csv'

    version = "1"
    nrandom = str(max(1000, random.random()*5000))
    xvar_psill = str(max(0.05,random.random()*1.0))
    minx = "-45.0"
    maxx = "45.0"
    miny = "-22.5"
    maxy = "22.5"

    var1_vrange = str(max(0.1, random.random()*5))
    var1_error = str(max(0.1, random.random()*5))
    prop_acc = str(max(0.1, random.random()*.9))
    var1_error_vrange = str(max(0.1, random.random()*5))
    mod_error_magnitude = str(max(0.05, random.random()*2))
    trt_prc = str(max(0.2, random.random()*0.5))
    theta = "1"
    beta = str(max(0.2, random.random()*5))
    spill_vrange = str(max(0.1, random.random()*5))
    spill_magnitude= str(max(0.05, random.random()*2))
    cal= str(max(0.25, random.random()*1))
    sample_size = str(max(0.1, random.random()*1))
    tree_split_lim=str(max(10, random.random()*50))
    mod_error_vrange=str(max(0.1, random.random()*5))
    xvar_error_psill = str(max(0.05,random.random()*1.0))
    mod_error_psill = str(max(0.05,random.random()*1.0))
    trt_spill_sill = str(max(0.05,random.random()*1.0))



    R_ret = check_output("Rscript " +
                        sim_path + " " +
                        version + " " +
                        nrandom+ " " +
                        xvar_psill+ " " +
                        minx+ " "+
                        maxx+ " "+
                        miny+ " "+
                        maxy+ " "+
                        var1_vrange+ " "+
                        var1_error+ " "+
                        prop_acc+ " "+
                        var1_error_vrange+ " "+
                        mod_error_magnitude+ " "+
                        trt_prc+ " "+
                        theta+ " "+
                        beta+ " "+
                        spill_vrange+ " "+
                        spill_magnitude+ " "+
                        cal+ " "+
                        sample_size+ " "+
                        tree_split_lim+ " "+
                        mod_error_vrange+ " "+
                        out_path+ " " +
                        xvar_error_psill+ " "+
                        mod_error_psill+ " "+
                        trt_spill_sill,
			stderr=STDOUT,
                        shell=True)

    print c
    print R_ret
    print "========================"

    c += size
