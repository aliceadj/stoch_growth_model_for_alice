function getPars(units, parametrization; abx = 0.0, ns=9.4810)

    NA = 602214085700000000519317.96307935

    if units == "uM"
        SFmolectoM_uM = 1e6 / (NA * 1e-15)
    elseif units == "molecs"
        SFmolectoM_uM = 1
    else  
        error("Please specify units as either 'uM' or 'molecs'")
    end 

    SF = 22000

    if parametrization == "PNAs"
        s0_val               = 1e4
        ns_val               = 0.5
        vt_val               = 726.0
        vm_val               = 5800.0
        Km_val               = 1.0e3 * SFmolectoM_uM 
        Kt_val               = 1.0e3 * SFmolectoM_uM 
        thetar_val           = (426.8693338968694) * SFmolectoM_uM 
        thetax_val           = (4.379733394834643) * SFmolectoM_uM 
        we_val               = 4.139172187824451 * SFmolectoM_uM 
        wr_val               = 929.9678874564831 * SFmolectoM_uM 
        wq_val               = 948.9349882947897 * SFmolectoM_uM   
        Kq_val               = 1.522190403737490e+05 * SFmolectoM_uM
        nq_val               = 4.
        gmax_val             = 1260.
        Kgamma_val           = 7 * SFmolectoM_uM #Kgamma
        nx_val               = 300.0
        nr_val               = 7549.0
        kb_val               = 1. / SFmolectoM_uM
        ku_val               = 1.0
        dm_val               = 0.1
        abx_val              = abx
        kon_val              = 0.00599
        Mref_val             = 1.0e8 * SFmolectoM_uM

  
    elseif parametrization == "NatComms"
        s0_val               = 1e6
        ns_val               = ns
        vt_val               = 660.0
        vm_val               = 5800.0
        Km_val               = 1.0e3 * SFmolectoM_uM 
        Kt_val               = 1.0e3 * SFmolectoM_uM 
        thetar_val           = (426.8693338968694*SF) * SFmolectoM_uM 
        thetax_val           = (4.379733394834643*SF) * SFmolectoM_uM 
        we_val               = 4.139172187824451 * 0.003 * SFmolectoM_uM 
        wr_val               = 929.9678874564831 * 0.008 * SFmolectoM_uM 
        wq_val               = 948.9349882947897 * 0.02 * SFmolectoM_uM 
        Kq_val               = 1.522190403737490e+05 * SFmolectoM_uM
        nq_val               = 4.
        gmax_val             = 12600.#1000
        Kgamma_val           = 7 * SF * SFmolectoM_uM #Kgamma
        nx_val               = 300.0
        nr_val               = 7549.0
        kb_val               = 1. / SFmolectoM_uM
        ku_val               = 1.0
        dm_val               = 0.1#1.
        abx_val              = abx 
        kon_val              = 0.0666#0.00599 # order of magnitude different? 
        Mref_val             = 1.0e8 * SFmolectoM_uM

    end 

    params = [:s0 => s0_val, :ns => ns_val, :vt => vt_val, :vm => vm_val, :Km => Km_val, :Kt => Kt_val, :thetar => thetar_val, :thetax => thetax_val,
               :we => we_val, :wr => wr_val, :wq => wq_val, :Kq => Kq_val, :nq => nq_val, :gmax => gmax_val, :Kgamma => Kgamma_val, :nx => nx_val, :nr => nr_val,
               :kb => kb_val, :ku => ku_val, :dm => dm_val, :abx => abx_val, :kon => kon_val, :M => Mref_val]


    return params
end

function getX0(model; ss=false, parametrization="NatComms")
    # NA = 602214085700000000519317.96307935
    # if units == "uM"
    #     SFmolectoM_uM = 1e6 / (NA * 1e-15)
    # elseif units == "molecs"
    #     SFmolectoM_uM = 1
    # else  
    #     error("Please specify units as either 'uM' or 'molecs'")
    # end 
    if ss
        if parametrization == "PNAs"
            init = [

                model.a => 12,#11.83,
                model.si => 143,#142.92,
                model.mm => 10,#9.92,
                model.mt => 10,#9.92,
                model.mq => 332,#331.63,
                model.mr => 36,#35.57,
                model.cm => 69,#69.01,
                model.ct => 69,#69.01,
                model.cq => 2307,#2306.51,
                model.cr => 802,#801.82,
                model.em => 7086,#7086.27,
                model.et => 7086,#7086.27,
                model.q => 2e5,#2.3683e5,
                model.r => 25,#25.4854
                model.zmm => 0.,
                model.zmt => 0.,
                model.zmq => 0.,
                model.zmr => 0.,
            ]
        elseif parametrization == "NatComms"
            init = [
                model.a => 3.168e8,#1.68e5,
                model.si => 676,#128,
                model.mm => 0,#0,
                model.mt => 0,#0,
                model.mq => 3,#1,
                model.mr => 1,#0,
                model.rmm => 1,#0,
                model.rmt => 1,#0,
                model.rmq => 433,#100,
                model.rmr => 544,#16,
                model.em => 1556,#2379,
                model.et => 2471,#2379,
                model.q => 507668,#2.8e5,
                model.r => 19388,#1709,
                model.zmm => 0,
                model.zmt => 0,
                model.zmq => 0,
                model.zmr => 0,
            ]
        else  
            error("Please specify units as either 'uM' or 'molecs'")
        end 

    else
        init = [

            model.a => 10., 
            model.si => 0.,
            model.mm => 0.,
            model.mt => 0.,
            model.mq => 0.,
            model.mr => 0.,
            model.cm => 0.,
            model.ct => 0.,
            model.cq => 0.,
            model.cr => 0.,
            model.em => 0.,
            model.et => 0.,
            model.q => 0.,
            model.r => 10.,
            model.zmm => 0.,
            model.zmt => 0.,
            model.zmq => 0.,
            model.zmr => 0.,
        ]
    end
    return init 
end 
