"""
Analysis of the results: plot data, and print tables
Note: Needs a set up PyCall environment
Tested with the follwowing python packages:
- tikzplotlib 0.9.8 
- numpy 1.20.2
- matplotlib 3.3.4
- seaborn 0.11.1
"""
#CHOOSE ONE OF THE FOLLOWING RUN TYPES
rtype = "ANAL_IMP"                      # produces performence plots from the electronic companion (Figure EC.1); approximation gaps, in-sample/out-of-sample stabilities (Figure 3); heuristic comparison (Figure 5) 
#rtype = "ANAL_OA_VS_BEN"                # produces Figure 4  
#rtype = "ANAL_EVALUATED_IMP"            # produces the 6 comparison plots (= main result) 
#rtype = "PRINT_TABLE_IMP"               # print tables from the electronic companion (Table EC.1 - EC.5)
#rtype = "PRINT_TABLE_CIMP"              # print table from the electronic companion (Table EC.6) 
#rtype = "ANAL_PROPAGATION_LENGTH_IMP"   # produces Figure 8 
#rtype = "ANAL_LEADER_UTILITY_CIMP"      # Figure EC.2 in electronic companion

# SET WORKING DIRECTORY
using Pkg
working_directory = dirname(@__FILE__)
cd(working_directory)

# ACTIVATE WORKING ENVIRONMENT (if not already active)
if Base.active_project() != dirname(@__FILE__) * "/Project.toml"
    pkg"activate ."
    pkg"instantiate"
end

# INCLUDES
include("declarations.jl")
include("plot_functions.jl")
include("helper_functions.jl")
include("print_tables.jl")
println("start script")


#set path and filename
rootpath = dirname(dirname(dirname(@__FILE__)))
path = resultpath = joinpath(rootpath, "results", "rawdata")
#imp
rtype == "ANAL_IMP"                     ? file = "20200822_IMP_RX.csv" : nothing
rtype == "ANAL_OA_VS_BEN"               ? file = "20210113_OA_BEN_w100.csv" : nothing
rtype == "ANAL_EVALUATED_IMP"           ? file = "20200904_IMP_diff_models_all_eval.csv" : nothing   #NOTE: this requires a dataframe with already evaluated objective functions
rtype == "ANAL_PROPAGATION_LENGTH_IMP"  ? file = "20200916_IMP_w100_all_models_wiener_prop_length.csv" : nothing
rtype == "PRINT_TABLE_IMP"              ? file = "20200916_IMP_w100_all_models.csv" : nothing
#cimp
rtype == "ANAL_LEADER_UTILITY_CIMP"     ? file = "20200916_CIMP_diff_leaderUtility.csv" : nothing
rtype == "PRINT_TABLE_CIMP"             ? file = "20200916_CIMP_diff_leaderUtility.csv" : nothing

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
if rtype ∉ ["ANAL_EVALUATED_IMP", "ANAL_EVALUATED_CIMP", "ANAL_PROPAGATION_LENGTH_IMP", "PRINT_TABLE_CIMP"]
    if rtype ∈ ["ANAL_IMP", "ANAL_OA_VS_BEN"]
        rdf.fractionOfResistantNodes = zeros( Float64, nrow(rdf) )                        
    end
    addEmptyDataColumns!( rdf )
    computeValues!(rdf, p)
    parseSeedSets!( rdf )
    parseSeedSetRanks!(rdf, :abs_rank_bc)
    computeSampleStability!( rdf )
end


#*****************************************
# ANALYIZE IMP
#*****************************************
if rtype == "ANAL_IMP"
    #PERFORMANCE PLOTS
    p.legendpos = "above"
    rename_instances!( rdf )

    #filter by scenarios
    df = rdf |> @filter((_.aggoutput == 0 ) .& (_.nsc == 100)) |> DataFrame
    plotRuntimeBySymbol( df, p, :ifile, "w=100" )
    df = rdf |> @filter((_.aggoutput == 0 ) .& (_.nsc == 250)) |> DataFrame
    plotRuntimeBySymbol( df, p, :ifile, "w=250" )
    df = rdf |> @filter((_.aggoutput == 0 ) .& (_.nsc == 500)) |> DataFrame
    plotRuntimeBySymbol( df, p, :ifile, "w=500" )
    df = rdf |> @filter((_.aggoutput == 0 ) .& (_.nsc == 750)) |> DataFrame
    plotRuntimeBySymbol( df, p, :ifile, "w=750" )

    
    # APPROXIMATION GAP / IN-SAMPLE /OUT-OF SAMPLE stabilities
    df = rdf |> @filter(_.aggoutput == 1 ) |> DataFrame
    xparams = [:nsc, :nsc, :nsc,
              ]
    yparams = [:gap_approx_rel, :inStability, :outStability,
              ]

    for i in eachindex( xparams )
        simpleSwarmBoxPlot( df, p, xparams[i], yparams[i], "R" )
    end


    # COMPARISON WITH HEURISTICS
    p.legendpos = "right"
    #simple
    df = rdf |> @filter((_.aggoutput == 1 ) .& (_.nsc == 100)) |> DataFrame
    stacked_df = getStackedDataFrameForComparisonOfObjectiveValuesWithHeursitic( df )
    xparams = [:methodname,
              ]
    yparams = [:obj_ratio,
              ]

    for i in eachindex( xparams )
        simpleSwarmPlot( stacked_df, p, xparams[i], yparams[i], "", :ifile )
    end

    df = DataFrame(stacked_df)
    df = df[:,[:ifile, :kF, :methodname, :obj_ratio, :ones]]
    for l in 1:nrow(df)
        df.methodname[l] == "MAR" ? df.methodname[l] = "MG" : nothing
    end
    cumulativePlotsBySymbols( df, p, :obj_ratio, :methodname )
end

if rtype == "ANAL_OA_VS_BEN"
    df = rdf |> @filter((_.aggoutput == 0 ) .& (_.nsc == 100)) |> DataFrame
    plotRuntimeBySymbol( df, p, :mtype )
end


if rtype in ["ANAL_PROPAGATION_LENGTH_IMP"]
    df_AGG = rdf |> @filter(_.aggoutput == 1) |> DataFrame

    xparams = [:mtypeString]
    yparams = [:propagation_length_mean_f]

    for i in eachindex(xparams)
        simpleBoxplot( df_AGG, p, xparams[i], yparams[i] )
    end
    df = df_AGG
end


if rtype in ["ANAL_EVALUATED_IMP"]
    parseSeedSetsAgain!(rdf)
    addEmptyDataColumnsAfterEvaluationOfObjectiveFunctions!( rdf )
    renameAndRoundDataFrameWithEvaluatedValues(rdf)
    reduceForwardingObjectiveValueBySeedSetSize!( rdf ) #NOTE: optional
    computeRelativeLossesWithRespectToDifferentModels!( rdf )
    parseNumberOfViewings!( rdf )
    computeMeanNumbersOfViewsAndForwards!(rdf)
    cleanOutlier!(rdf)
    computeSampleStability!( rdf )

    #plot gains only considering optimal solutions
    p.legendpos = "above"
    rename_instances!( rdf )
    opt_df          = removeNonOptimalSolutions(rdf)
    stacked_opt_df  = getStackedDataFrameForComparisonOfDifferentMethods( opt_df )
    my_methods      = unique(stacked_opt_df.mtypeString)
    my_kFs          = unique(stacked_opt_df.kF)
    my_kLs          = unique(stacked_opt_df.kL)
    sort!(my_methods)

    stacked_opt_df.value .= ((100 .- stacked_opt_df.value) ./ stacked_opt_df.value) .* 100

    for m in my_methods
        plot_df = stacked_opt_df |> @filter(_.mtypeString == m) |> DataFrame
        simpleSwarmBoxPlot( plot_df, p, :methodname, :value, m )
    end
end

if rtype in["PRINT_TABLE_IMP"]
    values_of_interest = [:ifile, :kF, :mtypeString, :rtCPLEX, :rtHbc, :rtHdc, :rtHed, :rtHMAR, :rtHpr, :rtHrm, :rtHtr, :obj_F_eval, :obj_Fbc_eval, :obj_Fdc_eval, :obj_Fed_eval, :obj_FMAR_eval, :obj_Fpr_eval, :obj_Frm_eval, :obj_Ftr_eval]
    rename_instances!( rdf )
    df_table = rdf |> @filter((_.aggoutput == 1) .& (_.nsc == 100))  |> DataFrame
    df_table = df_table[:, values_of_interest]
    orderdict = Dict(   "msg-college"       => 1,
                        "msg-email-eu"      => 2,
                        "soc-advogato"      => 3,
                        "soc-anybeat"       => 4,
                        "tw-austria"        => 5,
                        "tw-giftideas"      => 6,
                        "tw-greenenergy"    => 7,
                        "tw-naturelovers"   => 8,
                        "tw-organicfood"    => 9,
                        "tw-orms"           => 10,
                        "tw-skateboarding"  => 11,
                        "tw-travelling"     => 12,
                        5                   => 1,
                        10                  => 2,
                        15                  => 3,
                        "R50"               => 5,
                        "RX"                => 7,
                        "T"                 => 3,
                        "R25"               => 4,
                        "R75"               => 6,
                        "F"                 => 1,
                        "O"                 => 2 )
    #add old results with random resistance values
    file = "20200822_IMP_RX.csv"
    df_old = CSV.read( joinpath(path, file); copycols=true, delim=';' )
    df_old_AGG = df_old |> @filter((_.aggoutput == 1) .& (_.nsc == 100))  |> DataFrame
    df_old_AGG.mtypeString = Vector{String}(undef, nrow(df_old_AGG))
    for l=1:nrow(df_old_AGG)
        df_old_AGG.mtypeString[l] = "RX"
    end
    df_old_AGG = df_old_AGG[values_of_interest]
    rename_instances!( df_old_AGG )
    append!(df_table, df_old_AGG)

    #sort
    sort!(df_table, [:ifile, :kF, :mtypeString], by= x->orderdict[x])

    #round values
    values_to_round = [:rtCPLEX, :rtHbc, :rtHdc, :rtHed, :rtHMAR, :rtHpr, :rtHrm, :rtHtr, :obj_F_eval, :obj_Fbc_eval, :obj_Fdc_eval, :obj_Fed_eval, :obj_FMAR_eval, :obj_Fpr_eval, :obj_Frm_eval, :obj_Ftr_eval]
    roundValuesForTable!( df_table, values_to_round, 0 )
    df = df_table

    # Reduce objective value of method F by |F| for fair comparison.
    for l in 1:nrow(df)
        if df.mtypeString[l] == "F"
            df.obj_F_eval[l] -= df.kF[l] 
            df.obj_Fbc_eval[l] -= df.kF[l] 
            df.obj_Fdc_eval[l] -= df.kF[l] 
            df.obj_Fed_eval[l] -= df.kF[l] 
            df.obj_FMAR_eval[l] -= df.kF[l] 
            df.obj_Fpr_eval[l] -= df.kF[l] 
            df.obj_Frm_eval[l] -= df.kF[l] 
            df.obj_Ftr_eval[l] -= df.kF[l]
        end
    end

    #clean tuncRank and RM heuristic
    for l in 1:nrow(df)
        if df.ifile[l] == "msg-college"
            df.rtHrm[l]         = -969696.0
            df.obj_Frm_eval[l]  = -969696.0
            df.rtHtr[l]         = -969696.0
            df.obj_Ftr_eval[l]  = -969696.0
        end
        if df.ifile[l] == "msg-email-eu"
            df.rtHrm[l]         = -969696.0
            df.obj_Frm_eval[l]  = -969696.0
            df.rtHtr[l]         = -969696.0
            df.obj_Ftr_eval[l]  = -969696.0
        end
        if df.ifile[l] == "soc-advogato"
            df.rtHrm[l]         = -969696.0
            df.obj_Frm_eval[l]  = -969696.0
            df.rtHtr[l]         = -969696.0
            df.obj_Ftr_eval[l]  = -969696.0
        end
        if df.ifile[l] == "soc-anybeat"
            df.rtHrm[l]         = -969696.0
            df.obj_Frm_eval[l]  = -969696.0
            df.rtHtr[l]         = -969696.0
            df.obj_Ftr_eval[l]  = -969696.0
        end
    end

    #print table
    mydf = filter(row -> (row.ifile in ["tw-austria", "tw-giftideas"]), df)
    printLatexTableIMP(mydf, p)
    mydf = filter(row -> (row.ifile in ["tw-greenenergy", "tw-naturelovers", "tw-organicfood"]), df)
    printLatexTableIMP(mydf, p)
    mydf = filter(row -> (row.ifile in ["tw-orms", "tw-skateboarding", "tw-travelling"]), df)
    printLatexTableIMP(mydf, p)
    mydf = filter(row -> (row.ifile in ["msg-college", "msg-email-eu"]), df)
    printLatexTableIMP(mydf, p)
    mydf = filter(row -> (row.ifile in ["soc-advogato", "soc-anybeat"]), df)
    printLatexTableIMP(mydf, p)
end




#*****************************************
# ANALYIZE CIMP
#*****************************************
if rtype in ["ANAL_LEADER_UTILITY_CIMP"]
    convertLeaderUtilitiesToFactors( rdf )
    computeRelativeLossesWRTLeaderUtilityAndSetSimilarity( rdf )
    df_AGG = rdf |> @filter( (_.aggoutput == 1) .& (_.leaderUtilityFactor != 0) .& (_.fractionOfResistantNodes == 0.5)) |> DataFrame
    xparams = [:leaderUtilityFactor]
    yparams = [:loss_leaderUtility]
    for i in eachindex( xparams )
        simpleSwarmBoxPlot( df_AGG, p, xparams[i], yparams[i], "LU" )
    end
end

if rtype in["PRINT_TABLE_CIMP"]
    rename_instances!( rdf )
    values_of_interest = [:ifile, :kF, :leaderUtility, :rtCPLEX, :obj_F_eval]
    values_to_round    = [:rtCPLEX, :obj_F_eval]
    df_table = rdf |> @filter((_.aggoutput == 1) .& (_.nsc == 100))  |> DataFrame
    roundValuesForTable!( df_table, values_to_round, 0 )

    orderdict = Dict(   "msg-college"       => 9,
                        "msg-email-eu"      => 10,
                        "soc-advogato"      => 11,
                        "soc-anybeat"       => 12,
                        "tw-austria"        => 1,
                        "tw-giftideas"      => 2,
                        "tw-greenenergy"    => 3,
                        "tw-naturelovers"   => 4,
                        "tw-organicfood"    => 5,
                        "tw-orms"           => 6,
                        "tw-skateboarding"  => 7,
                        "tw-travelling"     => 8,
                        5                   => 1,
                        10                  => 2,
                        15                  => 3,
                        0                   => 1,
                        1                   => 2,
                        8                   => 3,
                        20                  => 4,
                        55                  => 5,
                        150                 => 6 )
    sort!(df_table, [:ifile, :kF, :leaderUtility], by= x->orderdict[x])
    df_table = df_table[:,values_of_interest]
    df_joined = getJoinedDataFrameWithRunTimesCIMP( df_table )
    df_joined = df_joined[:,[:ifile, :kF, :rtCPLEX, :rtCPLEX_1, :rtCPLEX_2, :rtCPLEX_3, :rtCPLEX_4, :rtCPLEX_5,
                                          :obj_F_eval, :obj_F_eval_1, :obj_F_eval_2, :obj_F_eval_3, :obj_F_eval_4, :obj_F_eval_5]]
    #df = df_table
    df = DataFrame(df_joined)

    #print table
    printLatexTableCIMP(df, p)
end

