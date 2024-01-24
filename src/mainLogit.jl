#******************************************************##
# MULTINOMIAL LOGIT (COMPETITIVE) INFLUENCE MAXIMIZATION
# Author: Michael Kahr
# Tested with: 
#   - Ubuntu 18.04 
#   - Julia 1.4, 
#   - CPLEX 12.9 and 12.10
#******************************************************

# PARAMETERS TO BE SPECIFIED BY USER
const CPLEX_DIR = "/opt/ibm/ILOG/CPLEX_Studio1210/cplex/bin/x86-64_linux"   # use either CPLEX 12.9, or 12.10

#SET WORKING DIRECTORY
using Pkg
working_directory = dirname(@__FILE__)
cd(working_directory)

# ACTIVATE WORKING ENVIRONMENT (if not already active)
if Base.active_project() != joinpath(dirname(@__FILE__), "Project.toml")
    ENV["CPLEX_STUDIO_BINARIES"] = CPLEX_DIR
    pkg"activate ."
    pkg"instantiate"
end

# INCLUDES & PACKAGES
using ArgParse
using LightGraphs, MetaGraphs
using Distributions
using Random
using TimerOutputs
using JuMP
using CPLEX
using MathOptInterface
const MOI = MathOptInterface
using StatsBase
using DataStructures
using SparseArrays
using Suppressor
include("var_declarations.jl")              # declaration of instance / results struct / and exception handler
include("params.jl")                        # parse program options
include("misc.jl")                          # miscellaneous functions
include("instance_reader.jl")               # instance reader
include("graph_functions.jl")               # several graph functions
include("heuristics.jl")                    # heuristics
include("cplex_stuff.jl")                   # cplex related stuff (e.g., parameters, getters, create empty model)
include("models_IMP.jl")                    # IMP models
include("models_CIMP.jl")                   # CIMP models
include("solution_checker.jl")              # check solution
include("write_output.jl")                  # write data to csv file


# MAIN
function main( params::Dict )
    println("Starting main...")
    
    # SET PATHS
    root_path = dirname(@__DIR__)
    params["root_path"] = root_path
    params["ipath"] = joinpath(root_path, "data", "instances")      # path to instances
    params["opath"] = joinpath(root_path, "results")                # path where results are written

    # DEFINE GLOBAL VARIABLES (used to compute average values after sample average approximation)
    #aggregated values (over all SSA iterations)
    aobj_val_F::Float64     = 0.0       #avg. objective value (follower)
    abest_bound::Float64    = 0.0       #avg. best bound
    anBB::Float64           = 0.0       #avg. number of BB nodes
    anLCUTS::Float64        = 0.0       #avg. number lazy cuts
    anUCUTS::Float64        = 0.0       #avg. number user cuts
    anBCUTS::Float64        = 0.0       #avg. number Benders cuts
    #aggregated runtimes (average over all SSA iterations)
    artCPLEX::Float64       = 0.0       #avg. runtime CPLEX
    artBM::Float64          = 0.0       #avg. runtime for building model
    artLG::Float64          = 0.0       #avg. runtime for creating live-graphs (i.e., scenario graphs)
    artR::Float64           = 0.0       #avg. runtime for creating activation sets
    artMNC::Float64         = 0.0       #avg. runtime for calculating set marginal node contributions
    artHMAR::Float64        = 0.0       #avg. runtime for marginal gain heuristic
    artCSOL::Float64        = 0.0       #avg. runtime for checking solution
    artLazyCB::Float64      = 0.0       #avg. runtime in lazy callback
    artUserCB::Float64      = 0.0       #avg. runtime in user callback
    artTOT::Float64         = 0.0       #avg. total runtime


    # GET / SET PARAMETERS ( cf., params.jl )
    # get params
    println("Parsed args:")
    for (arg,val) in params
        println( "$arg  =>  $val" )
    end
    if params["mtype"] in [1,2,3] #more precision for the linear case
        params["relcuttol"] = 1.001
    end
    isdir(params["opath"]) ? nothing : mkdir(params["opath"])

    # set random seed
    Random.seed!( params["rseed"] );

    # Create a TimerOutput object ( cf., pkg TimerOutputs )
    to::TimerOutput = TimerOutput()

    # true/false write output file
    wo::Bool = params["writeoutput"]

    # set test-instances if julia is not started externally (cf., misc.jl)
    params["console"] ? setInstancePathAndName!( params ) : nothing

    #set instance type
    if occursin( "D-", params["ifile"] )    # artificial testing instances     
        params["itype"] = 1
    end
    if occursin( ".im", params["ifile"] )   # from SNAP database
        params["itype"] = 2
    end
    if occursin( "TW-", params["ifile"] )   # twitter
        params["itype"] = 3
    end


    # CREATE INSTANCE / RESULTS / BFSdata  (cf., var_declarations.jl)
    inst::instance = createEmptyInstance( params )
    res::results   = createEmptyResults( params, inst )

    #sizehints:
    sizehint!( inst.SF, params["kF"] )
    sizehint!( inst.SL, params["kL"] )


    # READ INSTANCE ( cf., instance_reader.jl )
    res.stage = "read_instance"
    if( params["ifile"] == "D-n25-k4-b0.3-beta2.0-18.0-i1" )    # instance used for precompilation
        params["itype"] = 1
        wo = false
    end
    wo ? writeOutputToCSV( params, inst, res ) : nothing
    try
        if params["itype"] == 1     #artificial instances
            @timeit to "read_Instance" read_Artificial_Instance!( inst, params )
        elseif params["itype"] == 2 # snap instances
            @timeit to "read_Instance" read_Snap_Instance!( inst, params )
        elseif params["itype"] == 3 # twitter instances
            @timeit to "read_Instance" read_Twitter_Instance_MetaGraphs!( inst, params )
        end
        @assert( params["kF"] + params["kL"] <= inst.nnodes, "|SF|+|SL| > |V|")
    catch msg
        println( "####READ INSTANCE EXCEPTION: ", msg )
        wherefrom = "read instance"
        handleExceptions( inst, res, params, wo, wherefrom, msg )
    end
    inst.memlimit = params["memlimit"]

    #CREATE EMPTY BFS DATA
    bfs::bfsdata  = createEmptyBFSdata( inst.nnodes, inst.nsc )


    # COMPUTE HASH VALUE
    inst.hashval = hash( [inst.ifile, inst.mtype, inst.nsc, inst.kL, inst.kF] )

    # READ LEADERS SEED SET
    if inst.kL != 0
        read_SSL!( inst, params )

        # set logical seed nodes
        inst.inSL = falses( inst.nnodes )
        for i ∈ inst.SL
            inst.inSL[i] = true
        end
        print("done!\n")
    end

    # COMPUTE RESISTNACE VALUES (cf., graph_functions.jl)
    res.stage = "get_resistance_values"
    wo ? writeOutputToCSV( params, inst, res ) : nothing
    if inst.mtype in [1,2,3]
        params["fractionOfResistantNodes"] = 0.0
    end
    Random.seed!( 43 );
    getResistanceValues!( inst, params )


    #start sample average approximation (SAA)
    for mySSAItr = 1:inst.nSSAItr
        inst.memok = true
        Random.seed!( mySSAItr )
        res.mySSAItr = mySSAItr
        println("")
        println("******* SAA Iteeration $mySSAItr **************")

        # CREATE LIVE-GRAPHS (cf. graph_functions.jl )
        print("creating ", inst.nsc, " live-graphs...." )
        res.stage = "create-lg"
        wo ? writeOutputToCSV( params, inst, res ) : nothing
        try
            @timeit to "createLiveGraphs!" createLiveGraphs!( inst )
            if !inst.memok
                res.msg = "##mem too high creating livegraphs"
                @warn "##mem too high creating livegraphs,... continue SAA iteration"
                wo ? writeOutputToCSV( params, inst, res ) : nothing
                continue #SAA iteration
            end
            res.rtLG    = round( TimerOutputs.time(to["createLiveGraphs!"]) / 10^9 ; digits=2 )
            res.nLG     = getNElementsInLiveGraph( inst )
            if params["debug"] == 0
                res.memLG = ceil( Base.summarysize( inst.lgo ) / (1024 * 1024) ) + ceil( Base.summarysize( inst.lgi ) / (1024 * 1024) ) #memory consumption of variable (MiB)
            end
        catch msg
            println( "####CREATE LIVE-GRAPH EXCEPTION: ", msg )
            wherefrom = "create livegraph"
            handleExceptions( inst, res, params, wo, wherefrom, msg )
        end
        print("done! \n")


        # COMPUTE REVERSE ACTIVATION SETS (cf., graph_functions.jl)
        res.stage = "compute_T_R"
        wo ? writeOutputToCSV( params, inst, res ) : nothing
        print("compute A^-_ω(i)...")
        try
            if( inst.kL == 0 )  #IMP
                @timeit to "getValidActivators!" getValidActivatorsIMP!( inst, bfs )
            else                #CIMP
                if inst.mtype in [4,5]  # compute scenario dependent resistance values
                    @timeit to "getValidActivators!" compute_rw!( inst, bfs )
                end
                @timeit to "getValidActivators!" getValidActivatorsCIMP!( inst, bfs )
            end
            if !inst.memok
                res.msg = "##mem too high creating sets R"
                @warn "##mem too high creating sets R,... continue SAA iteration"
                wo ? writeOutputToCSV( params, inst, res ) : nothing
                continue #SAA iteration
            end
            res.nSingletons = sum( inst.isSingleton )
            res.rtR      = round( TimerOutputs.time(to["getValidActivators!"]) / 10^9 ; digits=2 )
            res.nR       = getNElementsInReachabilityset( inst )
            if params["debug"] == 0
                res.memR     = ceil( (Base.summarysize( inst.Rin ) + Base.summarysize( inst.Rout )) / (1024*1024) )
            end
        catch msg
            println( "####EXCEPTION computing A^-_omega(j) ", msg )
            wherefrom = "compute A^-_omega(j)"
            handleExceptions( inst, res, params, wo, wherefrom, msg )
        end
        print("done! \n")


        # START HEURISTIC (cf., heuristics.jl)
        res.stage = "heuristic"
        wo ? writeOutputToCSV( params, inst, res ) : nothing
        try
            if( params["htype"] == 3 )  #marginal gain
                print("start heuristic (type=3)...")
                if inst.kL == 0 #IMP
                    @timeit to "heurisitcMAR" inst.SFMAR  = heuristicMarginalGainIMP_all_models( inst, bfs )
                    res.obj_FMAR = computeObjectiveValuesFromSetIMP( inst, bfs, inst.SFMAR )             
                    res.rtHMAR  = round( TimerOutputs.time(to["heurisitcMAR"]) / 10^9 ; digits=2 )
                    push!(inst.SFMARSAA, inst.SFMAR)
                else
                    @timeit to "heurisitcMAR" inst.SFMAR  = heuristicMarginalGainCIMP_all_models( inst, bfs )
                    res.obj_FMAR = computeObjectiveValuesFromSetCIMP( inst, bfs, inst.SFMAR )       
                    res.rtHMAR  = round( TimerOutputs.time(to["heurisitcMAR"]) / 10^9 ; digits=2 )
                    push!(inst.SFMARSAA, inst.SFMAR)
                end
                if !inst.memok
                    res.msg = "##mem too high creating sets heristic MAR"
                    @warn "##mem too high creating sets heristic MAR"
                    wo ? writeOutputToCSV( params, inst, res ) : nothing
                end
                res.bestboundMARSAA[mySSAItr] = res.obj_FMAR
                print("done!\n")
                println("obj_FMAR: ", res.obj_FMAR)
            end
        catch msg
            println( "####SOLVE HEURISTIC EXCEPTION: ", msg )
            wherefrom = "heuristic"
            handleExceptions( inst, res, params, wo, wherefrom, msg )
        end

        #OUTDEGREE HEURISTIC
        expected_outdegree_heuristic!( inst, res )
        if inst.kL == 0
            res.obj_FDEG = computeObjectiveValuesFromSetIMP( inst, bfs, inst.SFDEG )
        else
            res.obj_FDEG = computeObjectiveValuesFromSetCIMP( inst, bfs, inst.SFDEG )
        end

        # ADJUST CPLEX TIMELIMIT
        timelimit::Int = 0
        timelimit = round(Int, params["timelimit"] - ( res.rtR ) )


        # CREATE CPLEX MODEL
        model = createEmptyModelMOI( params, timelimit )
        MOI.set(model, MOI.NumberOfThreads(), 1)


        # BUILD MODELS
        res.stage = "build_model"
        wo ? writeOutputToCSV( params, inst, res ) : nothing
        @assert( 0 < inst.mtype <=5, "inst.mtype out of bounds" )
        try
            if inst.kL == 0     #IMP
                if inst.mtype in [1,2,3,5]  #F, T, O, R
                    @timeit to "build_cpx-model" buildModel_BEN( inst, params, model, to )
                end
                if( inst.mtype == 4 ) #R with OA
                    @timeit to "build_cpx-model" buildModel_OA( inst, params, model )
                end
            else                #CIMP
                if inst.mtype in [1,2,3,5]   #F, T, O, R
                    @timeit to "build_cpx-model" buildModel_BENCIMP( inst, params, model )
                end
                if( inst.mtype == 4 ) #R with OA
                    @warn "Outer approximation not implemented for CIMP. Abort!"
                    return nothing
                end
            end

            res.rtBM = round( TimerOutputs.time(to["build_cpx-model"]) / 10^9 ; digits=2 )
        catch msg
            println( "####BUILD MODEL EXCEPTION: ", msg )
            wherefrom = "build model"
            handleExceptions( inst, res, params, wo, wherefrom, msg )
        end


        # CALLBACKS (cf., models.jl)
        # procstatus: used to abort solving if memory consumption exceeds set memory limit
        memorylimit::Int = params["memlimit"]
        function procStatusHeuristic( cb_data::CPLEX.CallbackContext )
            return_value = startProcStatus( cb_data, inst, memorylimit, res )                 # returns either CPLEX.terminate() or nothing
            return return_value
        end
        MOI.set(model, MOI.HeuristicCallback(), procStatusHeuristic)        #NOTE: only use for Linux

        # IMP
        # lazy GBD callback
        function lazyBendersOptimalityCutsIMP( cb_data::CPLEX.CallbackContext )
            @timeit to "LazyCB" startLazyBendersCallback( cb_data, inst, params, res, bfs, model )
            return nothing
        end

        # user GBD callback
        function userBendersOptimalityCutsIMP( cb_data::CPLEX.CallbackContext )
            @timeit to "UserCB" startUserBendersCallback( cb_data, inst, params, res, bfs, model, to )
            return nothing
        end

        # lazy OA callback
        function lazyOAOptimalityCutsIMP( cb_data::CPLEX.CallbackContext )
            wherefrom = "LCB"
            @timeit to "LazyCB" startOACallbackIMP( cb_data, inst, params, res, bfs, model, wherefrom )
            return nothing
        end

        # user OA callback
        function userOAOptimalityCutsIMP( cb_data::CPLEX.CallbackContext )
            wherefrom = "UCB"
            @timeit to "LazyCB" startOACallbackIMP( cb_data, inst, params, res, bfs, model, wherefrom ) 
            return nothing
        end

        # CIMP
        # lazy GBD callback
        function lazyBendersOptimalityCutsCIMP( cb_data::CPLEX.CallbackContext )
            @timeit to "LazyCB" startLazyBendersCallbackCIMP( cb_data, inst, params, res, bfs, model )
            return nothing
        end

        # user GBD callback
        function userBendersOptimalityCutsCIMP( cb_data::CPLEX.CallbackContext )
            @timeit to "UserCB" startUserBendersCallbackCIMP( cb_data, inst, params, res, bfs, model )
            return nothing
        end


        # add callbacks
        if inst.kL == 0 #IMP
            if params["mtype"] in [1,2,3,5]  # F, T, O, R
                MOI.set(model, MOI.UserCutCallback(), userBendersOptimalityCutsIMP)
                MOI.set(model, MOI.LazyConstraintCallback(), lazyBendersOptimalityCutsIMP)
            end
            if params["mtype"] == 4 # R with OA
                MOI.set(model, MOI.UserCutCallback(), userOAOptimalityCutsIMP)
                MOI.set(model, MOI.LazyConstraintCallback(), lazyOAOptimalityCutsIMP)
            end
        else            #CIMP
            if params["mtype"] in [1,2,3,5]  #F, T, O, R
                MOI.set(model, MOI.UserCutCallback(), userBendersOptimalityCutsCIMP)
                MOI.set(model, MOI.LazyConstraintCallback(), lazyBendersOptimalityCutsCIMP)
            end
        end


        # SOLVE
        print( "solve model...\n" )
        res.stage = "solve_model"
        res.bestbound = typemax( Float64 )
        wo ? writeOutputToCSV( params, inst, res ) : nothing

        # "typical" solving
        try
            @timeit to "solve(model)" optimize!( model )
            status = termination_status( model )
            println( "Status: ", status )
            res.cpxStatus = string( status )
        catch msg
            println( "####SOLVE MODEL EXCEPTION: ", msg )
            wherefrom = "solve model"
            res.cpxStatus = "NotSolved"
            res.exitflag = -1
            reset_timer!(to)
            handleExceptions( inst, res, params, wo, wherefrom, msg )
        end

        if result_count( model ) > 0
            #retrieve runtimes
            if (res.cpxStatus == "OPTIMAL") || (res.cpxStatus == "ALMOST_OPTIMAL") || (res.cpxStatus == "TIME_LIMIT") || (res.cpxStatus == "OTHER_ERROR")
                obj = round( objective_value( model ) ; digits=2 )
                res.rtCPLEX = round( TimerOutputs.time(to["solve(model)"]) / 10^9 ; digits=2 )
                try
                    res.rtLazyCB = round( TimerOutputs.time(to["solve(model)"]["LazyCB"]) / 10^9 ; digits=2 )
                catch
                    res.rtLazyCB = 0
                end
                try
                    res.rtUserCB = round( TimerOutputs.time(to["solve(model)"]["UserCB"]) / 10^9 ; digits=2 )
                catch
                    res.rtUserCB = 0
                end

                # get more values
                try
                    res.bestbound = round( MOI.get(model, MOI.ObjectiveBound()); digits=2 )
                    res.nBB       = MOI.get(model, MOI.NodeCount())
                    res.bestboundSAA[mySSAItr] = MOI.get(model, MOI.ObjectiveBound())
                    res.nLCUTSC   = cpx_get_num_cuts( model, CPLEX.CPX_CUT_TABLE )
                    res.nUCUTSC   = cpx_get_num_cuts( model, CPLEX.CPX_CUT_USER )
                    println("best bound: ", res.bestbound)
                catch msg
                    println( "####GET BEST BOUND EXCEPTION: ", msg )
                    wherefrom = "get best bound"
                    handleExceptions( inst, res, params, wo, wherefrom, msg )
                end
            end
        else
            #found no result
            res.rtCPLEX = params["timelimit"]
            res.cpxStatus = "NO_SOLUTION_FOUND"
        end


        # GET RESULTS & CHECK SOLUTION
        res.stage = "get_results"
        wo ? writeOutputToCSV( params, inst, res ) : nothing

        #cplex status / results
        if result_count( model ) > 0    #found at least some solution
            if (res.cpxStatus == "OPTIMAL") || (res.cpxStatus == "ALMOST_OPTIMAL" )      #:CPX_STAT_OPTIMAL, :CPXMIP_OPTIMAL, :CPXMIP_OPTIMAL_TOL
                res.obj_F   = objective_value( model )
                res.gap = round( MOI.get(model, MOI.RelativeGap()) * 100 ; digits=2 )               
                res.exitflag = 1
                @timeit to "check_solution" check_solution!( inst, res, bfs, model )
                res.rtCSOL = round( TimerOutputs.time(to["check_solution"]) / 10^9 ; digits=2 )
                println( "Objective value: ", res.obj_F )
                res.obj_F   = round( res.obj_F ; digits=2 )
                if(!params["solveRoot"])
                    println( "Gap: ", res.gap )
                end
            elseif (res.cpxStatus == "TIME_LIMIT") || (res.exitflag == 3)            #:TIME_LIM, :MIP_ABORT
                res.exitflag == 3 ? res.msg = "hit memlimit" : nothing
                if res.cpxStatus == "TIME_LIMIT"
                    res.exitflag = 2                                                   #hit TIME_LIM
                    res.msg = "hit timelimit"
                end

                wo ? writeOutputToCSV( params, inst, res ) : nothing                  
                if result_count( model ) == 0
                    res.nexceptions += 1
                else
                    try
                        res.obj_F     = round( objective_value( model ) ; digits=2 )
                        res.gap       = round( MOI.get(model, MOI.RelativeGap()) * 100 ; digits=2 )
                        @timeit to "check_solution" check_solution!( inst, res, bfs, model )
                        res.rtCSOL    = round( TimerOutputs.time(to["check_solution"]) / 10^9 ; digits=2 )
                        res.bestbound = round( MOI.get(model, MOI.ObjectiveBound()); digits=2 )
                        res.nBB       = MOI.get( model, MOI.NodeCount() )
                        println( "best bound: ", res.bestbound )
                        println( "Objective value: ", res.obj_F )
                        println( "Gap: ", res.gap )
                    catch msg
                        println( "#### ACCESS SOLUTION IN UserLimit EXCEPTION: ", msg )
                        println( "Objective value: NaN" )
                        println( "Gap: NaN" )
                        res.exitflag=6
                        res.msg = "access solution in UserLimit failed"
                        res.obj_F = -1
                        res.gap = -1
                        wo ? writeOutputToCSV( params, inst, res ) : nothing
                    end
                end
            else
                res.exitflag = 7                                                        # UNKNOWN STATE
                res.msg = res.cpxStatus
                @warn "non considered case in solving"
            end
        else #found no solution
            @warn "found no integer solution yet"
            res.exitflag=6
            res.msg = "found no integer solution yet"
            res.obj_F = -1
            res.gap = -1
        end
        if params["debug"] >= 2
            global tmodel = model
        end

        # WRITE OUTPUT ( cf., write_output.jl )
        if result_count( model ) > 0    #found at least some solution
            #COMPUTE LP RELAXATION GAP
            res.LPgap = ((res.LPrelax - res.obj_F_solcheck) / (0.0000000001 + res.obj_F_solcheck)) * 100

            #EVALUATE OBJECTIVE FUNCTION
            print("evaluating objective functions...")
            if !isempty(inst.SF)
                Random.seed!( 666 );
                res.obj_F_eval, res.obj_L_eval, u_vals = evaluateObjectiveFunction( inst, bfs, inst.SL, inst.SF )
                println("")
                println("obj_F_eval:", res.obj_F_eval)
                println("obj_L_eval:", res.obj_L_eval)
                Random.seed!( 666 );
                avg_min_f, avg_mean_f, avg_max_f = evaluatePropagationLength( inst, bfs, inst.SL, inst.SF, "F" )
                Random.seed!( 666 );
                avg_min_v, avg_mean_v, avg_max_v = evaluatePropagationLength( inst, bfs, inst.SL, inst.SF, "V" )
                println("f: min = $avg_min_f \tmean = $avg_mean_f \tmax = $avg_max_f")
                println("v: min = $avg_min_v \tmean = $avg_mean_v \tmax = $avg_max_v")
            end

            if !isempty(inst.SFMAR)
                Random.seed!( 666 );
                res.obj_FMAR_eval, res.obj_LMAR_eval, u_vals = evaluateObjectiveFunction( inst, bfs, inst.SL, inst.SFMAR )
                println("obj_FMAR : ", res.obj_FMAR)
                println("obj_FMAR_eval : ", res.obj_FMAR_eval)
            end
            println("done!")


            # WRITE OUTPUT ( cf., write_output.jl )
            res.rtTOT = round( TimerOutputs.tottime( to ) / 10^9 ; digits=2 )
            res.stage = "finished"
            wo ? writeOutputToCSV( params, inst, res ) : nothing

            #AGGREGATE DATA:
            if !isempty(inst.SF)
                push!(inst.SFSAA, inst.SF)
                aobj_val_F += res.obj_F
            else
                aobj_val_F += 0
            end

            abest_bound += res.bestbound
            anBB        += res.nBB
            anLCUTS     += res.nLCUTS
            anUCUTS     += res.nUCUTS
            anBCUTS     += anLCUTS + anUCUTS
            artCPLEX    += res.rtCPLEX
            artBM       += res.rtBM
            artLG       += res.rtLG
            artR        += res.rtR
            artMNC      += res.rtMNC
            artHMAR     += res.rtHMAR
            artCSOL     += res.rtCSOL
            artLazyCB   += res.rtLazyCB
            artUserCB   += res.rtUserCB
            artTOT      += res.rtTOT

            #CLEAR DATA
            print(to)
            if inst.nSSAItr > 1
                inst.SF = Int[]
                reset_timer!(to)
                reset_results!(res)
            end
            print("\n")
        end #if found solution

    end #SAAItr

    if inst.nSSAItr > 1
        # set param for aggregated file output
        res.aggoutput = 1
        res.mySSAItr  = 0

        #evaluate seed sets on large set of scenarios 
        println("Evaluating objective function:")
        best_set_id::Int = 1
        best_set_idIND::Int = 1
        best_obj_val::Float64 = 0.0
        best_obj_valL::Float64 = 0.0
        best_obj_valIND::Float64 = 0.0
        best_obj_valLIND::Float64 = 0.0
        dummy_counter::Int = 1
        best_u_vals = zeros(Float64, inst.Nsc)

        #evaluate solutions
        if !isempty( inst.SFSAA )
            for SF ∈ inst.SFSAA
                Random.seed!( 666 );
                @timeit to "eval" obj_F_eval, obj_L_eval, u_vals = evaluateObjectiveFunction( inst, bfs, inst.SL, SF )
                rt = round( TimerOutputs.time(to["eval"]) / 10^9 ; digits=2 )
                println("eval SS$dummy_counter: $rt s")
                if obj_F_eval > best_obj_val
                    best_set_id   = dummy_counter
                    best_obj_val  = obj_F_eval
                    best_obj_valL = obj_L_eval
                    best_u_vals   = deepcopy(u_vals)
                end
                dummy_counter += 1
            end
            res.rtEVAL = round( TimerOutputs.time(to["eval"]) / 10^9 ; digits=2 )
            println("found best solution at pos: $best_set_id")
            res.best_ss_id = best_set_id
            inst.SF = inst.SFSAA[best_set_id]
            res.obj_F_eval = best_obj_val
            res.obj_L_eval = best_obj_valL
            res.gap_approx_rel = computeApproximationGap( inst, res, best_u_vals, res.bestboundSAA )
            if(params["debug"] >= 2)
                global tbest_u_vals = best_u_vals
            end
        end

        #evaluate MAR heuristic solutions
        if !isempty(inst.SFMARSAA)
            dummy_counter = 1
            for SF ∈ inst.SFMARSAA
                Random.seed!( 666 );
                @timeit to "evalIND" obj_FMAR_eval, obj_LMAR_eval, u_vals = evaluateObjectiveFunction( inst, bfs, inst.SL, SF )
                rt = round( TimerOutputs.time(to["evalIND"]) / 10^9 ; digits=2 )
                println("eval SS$dummy_counter: $rt s")
                if obj_FMAR_eval > best_obj_valIND
                    best_set_idIND   = dummy_counter
                    best_obj_valIND  = obj_FMAR_eval
                    best_obj_valLIND = obj_LMAR_eval
                    best_u_vals      = deepcopy(u_vals)
                end
                dummy_counter += 1
            end
            res.rtEVALMAR = round( TimerOutputs.time(to["evalIND"]) / 10^9 ; digits=2 )
            res.best_ss_idMAR = best_set_idIND
            inst.SFMAR = inst.SFMARSAA[best_set_idIND]
            res.obj_FMAR_eval = best_obj_valIND
            res.obj_LMAR_eval = best_obj_valLIND
            println("found best solution MAR at pos: $best_set_idIND")
            res.gap_approx_relMAR = computeApproximationGap( inst, res, best_u_vals, res.bestboundMARSAA )
        end


        #COMPUTE AVERAGE VALUES
        res.obj_F       = aobj_val_F  / length(inst.SFSAA)
        res.bestbound   = abest_bound / (inst.nSSAItr - res.nexceptions)
        res.nBB         = round(Int, anBB        / (inst.nSSAItr - res.nexceptions))
        res.nLCUTS      = round(Int, anLCUTS     / (inst.nSSAItr - res.nexceptions))
        res.nUCUTS      = round(Int, anUCUTS     / (inst.nSSAItr - res.nexceptions))
        res.nBENCUTS    = round(Int, anBCUTS     / (inst.nSSAItr - res.nexceptions))
        res.rtCPLEX     = artCPLEX   / (inst.nSSAItr - res.nexceptions)
        res.rtBM        = artBM      / (inst.nSSAItr - res.nexceptions)
        res.rtLG        = artLG      / (inst.nSSAItr - res.nexceptions)
        res.rtR         = artR       / (inst.nSSAItr - res.nexceptions)
        res.rtMNC       = artMNC     / (inst.nSSAItr - res.nexceptions)
        res.rtCSOL      = artCSOL    / (inst.nSSAItr - res.nexceptions)
        res.rtLazyCB    = artLazyCB  / (inst.nSSAItr - res.nexceptions)
        res.rtUserCB    = artUserCB  / (inst.nSSAItr - res.nexceptions)
        res.rtTOT       = artTOT     / (inst.nSSAItr - res.nexceptions)

    end

    #COMPUTE RANKINGS
    if( params["computeRankings"] )
        @timeit to "computeRankings" begin
            if !isempty(inst.SF)
                computeRankings!( inst, res, inst.SF, "BEN", params, to )
            end
            if !isempty(inst.SFMAR)
                computeRankings!( inst, res, inst.SFMAR, "MAR", params, to )
            end
        end
        Random.seed!( 666 );
        res.obj_Fbc_eval, res.obj_Lbc_eval, u_vals = evaluateObjectiveFunction( inst, bfs, inst.SL, inst.SFbc )
        Random.seed!( 666 );
        res.obj_Fdc_eval, res.obj_Ldc_eval, u_vals = evaluateObjectiveFunction( inst, bfs, inst.SL, inst.SFdc )
        Random.seed!( 666 );
        res.obj_Fcc_eval, res.obj_Lcc_eval, u_vals = evaluateObjectiveFunction( inst, bfs, inst.SL, inst.SFcc )
        Random.seed!( 666 );
        res.obj_Fed_eval, res.obj_Led_eval, u_vals = evaluateObjectiveFunction( inst, bfs, inst.SL, inst.SFed )
        Random.seed!( 666 );
        res.obj_Fpr_eval, res.obj_Lpr_eval, u_vals = evaluateObjectiveFunction( inst, bfs, inst.SL, inst.SFpr )
        Random.seed!( 666 );
        res.obj_Ftr_eval, res.obj_Ltr_eval, u_vals = evaluateObjectiveFunction( inst, bfs, inst.SL, inst.SFtr )
        if params["itype"] == 3
            Random.seed!( 666 );
            res.obj_Frm_eval, res.obj_Lrm_eval, u_vals = evaluateObjectiveFunction( inst, bfs, inst.SL, inst.SFrm )
        end
    end

    if isempty(inst.SF)
        res.msg = "no solution found"
    end

    # FINALIZE ( cf., write_output.jl )
    res.rtTOT = round( TimerOutputs.tottime( to ) / 10^9 ; digits=2 )
    res.stage = "finished"
    wo ? writeOutputToCSV( params, inst, res ) : nothing
    if( params["computeSL"] )
        wo ? writeOutputLeadersSeedSet( params, inst, res ) : nothing
        Random.seed!( params["rseed"] )
    end

    # PRINT RUNTIME/MEMORY STATISTICS
    println("")
    println("LPrelax: ", res.LPrelax )
    println("obj_L: ", res.obj_L )
    println("nBENCUTS: ", res.nBENCUTS )
    println("nLCUTS: ", res.nLCUTS )
    println("nUCUTS: ", res.nUCUTS )
    println("nBB: ", res.nBB )


    # give access to structs in REPL
    println("")
    if(params["debug"] >= 2)                            
        global tinst = inst
#        global tmodel = model
        global tres = res
        global tto = to
        global tbfs = bfs
        global tparams = params
    end

    println("exit program")

    return nothing
end #main()


function startMyMain()
    params = parse_commandline()
#    println("Parsed args:")
    for (arg,val) in params
#        println( "$arg  =>  $val" )
    end

    # start on local console (repl) / or on cluster
    if params["console"] == true
        setInstancePathAndName!( params )
        main( params )
        return
    end

    if params["precompile"] == false    #start main with args
        main( params )
    else                                #run precompile
        #set dummy params
        pre_params = deepcopy( params )
        pre_params["ifile"]         = "D-n25-k4-b0.3-beta2.0-18.0-i1"
        pre_params["itype"]         = 1
        pre_params["mtype"]         = 5
        pre_params["nscenarios"]    = 10
        pre_params["Nscenarios"]    = 10
        pre_params["nSSAItr"]       = 1
        pre_params["benf"]          = true
        pre_params["kL"]            = 0
        pre_params["kF"]            = 2
        pre_params["relcuttol"]     = 1.01
        pre_params["relcuttolU"]    = 1.01
        pre_params["prepro"]        = true
        pre_params["debug"]         = 0
        pre_params["cutType"]       = 5
        pre_params["fractionOfResistantNodes"] = 0.1
        pre_params["probV"]         = 0.5
        pre_params["solveRoot"]     = false
        pre_params["rootCuts"]      = false
        pre_params["rootCutsAlpha"] = 0.5
        pre_params["rootCutsLambda"]= 0.5
        pre_params["leaderUtility"] = 0
        pre_params["ssltype"]       = 4
        pre_params["computeSL"]     = true
        pre_params["timelimit"]     = 100
        pre_params["htype"]         = 0
        pre_params["hcallback"]     = false
        pre_params["OAtype"]        = 0
        pre_params["LPMethod"]      = 0

        println("\n##PRECOMPILE RUN...")
        @suppress begin
            main(pre_params)
        end

        println("\n##HOT RUN...")
        main(params)
    end

    return
end
startMyMain()
