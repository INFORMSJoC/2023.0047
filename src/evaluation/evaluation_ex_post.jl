"""
Ex post evaluation of seed sets on different metrics (used in Section 7.5)
Ex post evaluation of the average propagation length
"""
#choose one of the following run types
rtype = "EVALUATE_IMP_MODELS"                # cross-evaluation of seed sets on different metrics 
#rtype = "EVALUATE_CIMP_MODELS"              # produces Figure 4  
#rtype = "EVALUATE_PROPAGATION_LENGTH_IMP"   # evaluation of average propagation length

using Pkg
rootpath = dirname(dirname(@__DIR__))
working_directory = dirname(@__FILE__)
cd(working_directory)

# ACTIVATE WORKING ENVIRONMENT (if not already active)
if Base.active_project() != dirname(@__FILE__) * "/Project.toml"
    pkg"activate ."
    pkg"instantiate"
end

using LightGraphs, MetaGraphs
using Random
using TimerOutputs
include("declarations.jl")
include("helper_functions.jl")
include("../var_declarations.jl")
include("../instance_reader.jl")
include("../graph_functions.jl")
include("../params.jl")
println("start script")

#set path and filename
rootpath = dirname(dirname(dirname(@__FILE__)))
path = resultpath = joinpath(rootpath, "results", "raw_data")

# LOAD RAW DATA (mutable)
rdf = CSV.read( joinpath(resultpath, file); copycols=true, delim=';' )     #raw data

# SET UP PARAMETERS
p = parameters( 7200,
                #rdf[:tlim][1],
                opath,
                size( rdf, 1 ),
                unique(rdf[!,:ifile]),
                colvec,
                [[1,1,1,1], [10,0,10,0], [1,3,2,3], [5,2,0,2], [4,2,4,2], [6,1,6,1], [2,1,4,1], [4,1,3,1], [2,1,2,1], [4,0,5,1], [1,0,1,1], [1,2,3,4]],
                1.05,       #line thickness
                mydict,     #converts a column name to a string
                "serif",    #font type of figure labels
                5,          #granularity performance plot
                0.01,       #granularity gap plots
                false,       #save tikz
                true,       #show fliers
                "linear",   #performace plot type options {"linear", "log", "symlog", "logit", ...}
                "auto",     #position of the legend
                true,       #dodge in swarmplot
                0,          #xlim LB
                100,        #xlim UB
                0,          #ylim LB
                100 )       #ylim UB


#*****************************************
# PRECOMPUTATION STEPS
#*****************************************                
if rtype ∉ ["EVALUATE_PROPAGATION_LENGTH_IMP"]
    if rtype ∈ ["ANAL_IMP", "ANAL_OA_VS_BEN"]
        rdf.fractionOfResistantNodes = zeros( Float64, nrow(rdf) )                        
    end
    addEmptyDataColumns!( rdf )
    computeValues!(rdf, p)
    parseSeedSets!( rdf )
    parseSeedSetRanks!(rdf, :abs_rank_bc)
    computeSetSimilarityOverSAAIterations!( rdf )
    computeSampleStability!( rdf )
end

#*****************************************
# EX POST EVALUATION
#*****************************************    
#evaluate seed sets from different models
function evaluateSeedSets!( rdf )
    df_AGG = rdf |> @filter(_.aggoutput == 1 ) |> DataFrame
    fr_vec = unique(df_AGG.fractionOfResistantNodes)
    params = parse_commandline( ARGS )
    params["--Nscenarios"] = 100000

    # group by ifile
    for subdf_file in groupby(df_AGG, :ifile)
        #read instance
        params["ifile"] = subdf_file.ifile[1]
        if params["ifile"] in ["msg-college.im", "msg-email-eu.im", "soc-advogato.im", "soc-anybeat.im"]
            params["itype"] = 2
            params["ipath"] = joinpath(root_path, "data", "instances")
        else
            params["itype"] = 3
            params["ipath"] = joinpath(root_path, "data", "instances")
        end

        Random.seed!( params["rseed"] );
        inst::instance  = createEmptyInstance( params )
        inst.nsc        = 100
        inst.Nsc        = 100000
        params["itype"] == 2 ? read_Snap_Instance!( inst, params ) : nothing
        params["itype"] == 3 ? read_Twitter_Instance_MetaGraphs!( inst, params ) : nothing
        bfs::bfsdata    = createEmptyBFSdata( inst.nnodes, inst.nsc )

        #group by kL
        for subdf_kL in groupby(subdf_file, :kL)
            inst.kL     = subdf_kL.kL[1]
            params["kL"]= inst.kL

            #groupby kF
            for subdf_kF in groupby(subdf_kL, :kF)
                inst.kF     = subdf_kF.kF[1]
                params["kF"]= inst.kF

                #group by mtype
                for subdf_mtype in groupby(subdf_kF, :mtype)
                    inst.mtype      = subdf_mtype.mtype[1]
                    params["mtype"] = inst.mtype

                    if inst.mtype == 5
                        #group by fraction of resistant nodes
                        for subdf_fractionOfResistantNodes in groupby(subdf_mtype, :fractionOfResistantNodes)
                            if nrow(subdf_fractionOfResistantNodes) > 1
                                @warn "subdf_fraction has more than one row!!!"
                                #global foo = deepcopy()
                            end
                            #println("rowID: ", subdf_fractionOfResistantNodes.rowID[1])
                            rowID = subdf_fractionOfResistantNodes.rowID[1]
                            params["fractionOfResistantNodes"] = subdf_fractionOfResistantNodes.fractionOfResistantNodes[1]
                            #print(subdf_fractionOfResistantNodes)
                            println("evaluateing...\t", params["ifile"],"\tm=", params["mtype"], "\tkL=", params["kL"], "\tkF=", params["kF"], "\tfr=", params["fractionOfResistantNodes"])
                            SF = subdf_fractionOfResistantNodes.SFv[1]
                            SL = subdf_fractionOfResistantNodes.SLv[1]
                            for fr in fr_vec
                                params["fractionOfResistantNodes"] = fr
                                Random.seed!( 43 );
                                getResistanceValues!( inst, params )

                                Random.seed!( 666 );
                                uF_resistant, uL_resistant, uF_total, uL_total, uF_organic, uL_organic, uF_forwards, uL_forwards = evaluateObjectiveFunctionExPost( inst, bfs, SL, SF )
                                if fr == 0.25
                                    rdf.obj_resistant_F25[rowID] = uF_resistant
                                    rdf.obj_resistant_L25[rowID] = uL_resistant
                                end
                                if fr == 0.5
                                    rdf.obj_resistant_F50[rowID] = uF_resistant
                                    rdf.obj_resistant_L50[rowID] = uL_resistant
                                end
                                if fr == 0.75
                                    rdf.obj_resistant_F75[rowID] = uF_resistant
                                    rdf.obj_resistant_L75[rowID] = uL_resistant
                                end
                                rdf.obj_total_F[rowID]      = uF_total
                                rdf.obj_total_L[rowID]      = uL_total
                                rdf.obj_organic_F[rowID]    = uF_organic
                                rdf.obj_organic_L[rowID]    = uL_organic
                                rdf.obj_forwards_F[rowID]   = uF_forwards
                                rdf.obj_forwards_L[rowID]   = uL_forwards
                                rdf.nTotalViewsF[rowID]     = inst.nTotalViewsF
                                rdf.nOrganicViewsF[rowID]   = inst.nOrganicViewsF
                                rdf.nForwardsF[rowID]       = inst.nForwardsF
                                rdf.nTotalViewsL[rowID]     = inst.nTotalViewsL
                                rdf.nOrganicViewsL[rowID]   = inst.nOrganicViewsL
                                rdf.nForwardsL[rowID]       = inst.nForwardsL
                            end
                        end
                    else #mtype ∈ [1,2,3]
                        if nrow(subdf_mtype) > 1
                            @warn "subdf_mtype has more than one row!!!"
                            #global fut = subdf_mtype
                        end
                        #println("rowID: ", subdf_mtype.rowID[1])
                        rowID = subdf_mtype.rowID[1]
                        params["fractionOfResistantNodes"] = 0.0
                        println("evaluateing...\t", params["ifile"],"\tm=", params["mtype"], "\tkL=", params["kL"], "\tkF=", params["kF"], "\tfr=", params["fractionOfResistantNodes"])
                        SF = subdf_mtype.SFv[1]
                        SL = subdf_mtype.SLv[1]
                        #println("SF: ", SF)
                        #println("SL: ", SL)
                        for fr in fr_vec
                            params["fractionOfResistantNodes"] = fr
                            Random.seed!( 43 );
                            getResistanceValues!( inst, params )

                            Random.seed!( 666 );
                            uF_resistant, uL_resistant, uF_total, uL_total, uF_organic, uL_organic, uF_forwards, uL_forwards = evaluateObjectiveFunctionExPost( inst, bfs, SL, SF )
                            if fr == 0.25
                                rdf.obj_resistant_F25[rowID] = uF_resistant
                                rdf.obj_resistant_L25[rowID] = uL_resistant
                            end
                            if fr == 0.5
                                rdf.obj_resistant_F50[rowID] = uF_resistant
                                rdf.obj_resistant_L50[rowID] = uL_resistant
                            end
                            if fr == 0.75
                                rdf.obj_resistant_F75[rowID] = uF_resistant
                                rdf.obj_resistant_L75[rowID] = uL_resistant
                            end
                            rdf.obj_total_F[rowID]      = uF_total
                            rdf.obj_total_L[rowID]      = uL_total
                            rdf.obj_organic_F[rowID]    = uF_organic
                            rdf.obj_organic_L[rowID]    = uL_organic
                            rdf.obj_forwards_F[rowID]   = uF_forwards
                            rdf.obj_forwards_L[rowID]   = uL_forwards
                            rdf.nTotalViewsF[rowID]     = inst.nTotalViewsF
                            rdf.nOrganicViewsF[rowID]   = inst.nOrganicViewsF
                            rdf.nForwardsF[rowID]       = inst.nForwardsF
                            rdf.nTotalViewsL[rowID]     = inst.nTotalViewsL
                            rdf.nOrganicViewsL[rowID]   = inst.nOrganicViewsL
                            rdf.nForwardsL[rowID]       = inst.nForwardsL
                        end
                    end
                end
            end
        end
    end

    return
end


#evaluate seed sets from different models
function evaluatePropagationLength!( rdf )
    df_AGG = rdf |> @filter(_.aggoutput == 1 ) |> DataFrame
#    df_AGG = rdf |> @filter( (_.aggoutput == 1) .& (_.ifile == "msg-college.im") .& (_.kF == 5) ) |> DataFrame
#    df_AGG = rdf |> @filter( (_.aggoutput == 1) .& (_.ifile == "TW-skateboarding-20200810-Y2020") .& (_.kF == 5) ) |> DataFrame

    params = parse_commandline( ARGS )
    params["--Nscenarios"] = 100

    for l in 1:nrow(df_AGG)
        #read instance
        params["ifile"] = df_AGG.ifile[l]
        if params["ifile"] in ["msg-college.im", "msg-email-eu.im", "soc-advogato.im", "soc-anybeat.im"]
            params["itype"] = 2
            params["ipath"] = joinpath(root_path, "data", "instances")
        else
            params["itype"] = 3
            params["ipath"] = joinpath(root_path, "data", "instances")
        end

        Random.seed!( params["rseed"] );
        inst::instance  = createEmptyInstance( params )
        inst.nsc        = 100
        inst.Nsc        = 10000
        params["itype"] == 2 ? read_Snap_Instance!( inst, params ) : nothing
        params["itype"] == 3 ? read_Twitter_Instance_MetaGraphs!( inst, params ) : nothing
        bfs::bfsdata    = createEmptyBFSdata( inst.nnodes, inst.nsc )

        #set parameters
        inst.kL     = params["kL"]      = df_AGG.kL[l]
        inst.kF     = params["kF"]      = df_AGG.kF[l]
        inst.mtype  = params["mtype"]   = df_AGG.mtype[l]
        inst.SL     = df_AGG.SLv[l]
        inst.SF     = df_AGG.SFv[l]
        rowID       = df_AGG.rowID[l]

        println("evaluateing...\t", params["ifile"],"\tm=", df_AGG.mtypeString[l], "\tkL=", params["kL"], "\tkF=", params["kF"])

        #propagation length
        Random.seed!( 666 );
        @time rdf.propagation_length_min_f[rowID], rdf.propagation_length_mean_f[rowID], rdf.propagation_length_max_f[rowID] = evaluatePropagationLength( inst, bfs, inst.SL, inst.SF, "F" )
        println("f:, min = ", rdf.propagation_length_min_f[rowID], "\tmean = ", rdf.propagation_length_mean_f[rowID], "\tmax = ", rdf.propagation_length_max_f[rowID])
        Random.seed!( 666 );
        @time rdf.propagation_length_min_v[rowID], rdf.propagation_length_mean_v[rowID], rdf.propagation_length_max_v[rowID] = evaluatePropagationLength( inst, bfs, inst.SL, inst.SF, "V" )
        println("f:, min = ", rdf.propagation_length_min_v[rowID], "\tmean = ", rdf.propagation_length_mean_v[rowID], "\tmax = ", rdf.propagation_length_max_v[rowID])

        #write output
        mypath = resultpath
        myfile = "20200916_IMP_w100_all_models_wiener.csv"
#        CSV.write( joinpath(mypath, myfile * "_prop_length.csv"), rdf;  delim=';' )
    end

    return
end


if rtype in ["EVALUATE_IMP_MODELS", "EVALUATE_CIMP_MODELS"]
    computeSeedSetSimilarityOverDifferentModels!( rdf )
    evaluateSeedSets!( rdf )
    roundValuesOfInterest!( rdf )
    print("writing csv...")
#    CSV.write( path * file * "_eval.csv", rdf;  delim=';' )
    print("done!\n")
end


if rtype in ["EVALUATE_PROPAGATION_LENGTH_IMP"]
    rdf.propagation_length_min_f          = zeros( Float64, nrow(rdf) )
    rdf.propagation_length_mean_f         = zeros( Float64, nrow(rdf) )
    rdf.propagation_length_max_f          = zeros( Float64, nrow(rdf) )
    rdf.propagation_length_min_v          = zeros( Float64, nrow(rdf) )
    rdf.propagation_length_mean_v         = zeros( Float64, nrow(rdf) )
    rdf.propagation_length_max_v          = zeros( Float64, nrow(rdf) )
    rdf.wiener_index_v          = zeros( Float64, nrow(rdf) )                         #Wiener index viewing
    rdf.SFv                     = Vector{Vector{Int}}(undef, nrow(rdf))               #Followers Seed Set (converted from String to Vector)
    rdf.SLv                     = Vector{Vector{Int}}(undef, nrow(rdf))               #Followers Seed Set (converted from String to Vector)
    for l in nrow(rdf)
        rdf.SFv[l] = Int[]
        rdf.SLv[l] = Int[]
    end
    parseSeedSets!( rdf )
    evaluatePropagationLength!(rdf)  
end
