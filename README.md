### TNC-PilotProject
Repository for UVic/Foundry Spatial/TNC pilot project in California.

## Some notes on the various scripts:
There is a rough (but not prescribed) order in which the various models
contained in this directory should be run because initial conditions and
model parameters are 'daisy-chained' from one to the next. For instance, 
the steady-state output is used as initial conditions and to provide parameter
values for the transient model.

Here is an overview of the recommended sequence and relevant scripts:
1. Steady-state model
    - MODFLOW_Navarro-SteadyState.py = create and run steady-state model
    - MODFLOW_HTC_Navarro_SetUpSteadyStateWithPumping.py = create subdirectories 
    needed to run pumping scenarios on a high-throughput computing system. MODFLOW_Navarro-SteadyState.py
    must be run *on your HTC machine* before running this script.
        - launch_allRuns.sl = script to deploy all pumping runs (`sbatch launch_allRuns.sl`)
        - launch_thisRun_(version).sh = script to launch an individual model run, 
        which is copied into each subdirectory
        - Navarro-SteadyState_Template_(version).nam = template for NAM file so 
        that subdirectories can share input files.
    - MODFLOW_HTC_CheckFailure.py = check whether any models failed to converge.
    - MODFLOW_HTC_Navarro_(stream BC)-SummarizeLeakage.py = postprocess and inspect 
    MODFLOW model output.
    - MODFLOW_Navarro-SteadyState-WithPumping.py = pump an individual well and look 
    at the output. MODFLOW_Navarro-SteadyState.py must be run on the machine before
    running this script.
2. Transient model
    - MODFLOW_Navarro-Transient-SpinUp.py = run a long spin-up simulation. MODFLOW_Navarro-SteadyState.py
    must be run on the machine before running this script.
    - MODFLOW_HTC_Navarro_SetUpTransientWithPumping.py = create subdirectories needed to run
    pumping scenarios on a HTC system. MODFLOW_Navarro-Transient-SpinUp.py must be run before
    running this script (can be run on a different machine, as the appropriate outputs are tracked).
        - launch_thisRun.sh = script to launch an individual model run, which is copied into each subdirectory
        - postprocess_thisRun_(stream BC).py = script to postprocess an individual model run, 
        which is copied into each subdirectory
        - launch_allRuns.sl = script to deploy all pumping runs (`sbatch launch_allRuns.sl`), which
        includes postprocessing
        - postprocess_allRuns.sl = script to deploy all postprocessing; this shouldn't be necessary
        since launch_allRuns.sl now includes postprocessing.
    - MODFLOW_HTC_Navarro_SetUpTransientWithPumping-AddMoreRuns.py = supplement the subdirectories
    created with MODFLOW_HTC_Navarro_SetUpTransientWithPumping.py if there are additional wells you
    want to create and run.
    - MODFLOW_HTC_CheckFailure.py = check whether any models failed to converge.
    - MODFLOW_HTC_Navarro_CalculateCaptureFraction.R = Calculate the capture fraction based on the 
    postprocessed output from the HTC runs.
3. Analytical model
    - Navarro_DepletionApportionment.R = calculate depletion apportionment equations
        - \_MaskDryStreams = apply a mask of dry streams to the steady-state model to only apportion
        depletion into flowing reaches. This is only relevant when (i) using the SFR package AND
        (ii) not including any overland flow.
    - Navarro_Analytical_Transient.R = calculate streamflow depletion for each stream reach and pumping
    well. Needs output from Navarro_DepletionApportionment.R
        - \_DetermineSimulationLength.R = estimate analytical depletion for a long time to figure out 
        when it asymptotes. Needs output from Navarro_DepletionApportionment.R
4. Model comparison
    - Navarro_CompareMODFLOWtoDepletionApportionment_SteadyState.R = compare MODFLOW and analytical results
    for steady-state model.
    - Navarro_CompareMODFLOWtoDepletionApportionment_Transient.R = compare MODFLOW and analytical results
    for transient model.

        
    
    
