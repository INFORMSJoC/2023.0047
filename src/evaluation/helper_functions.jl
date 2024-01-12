
# ADD NEW COLS
function addEmptyDataColumns!( rdf::DataFrame )
    rdf.rt                      = zeros( Float64, nrow(rdf) )                         #aggregated runtime
    rdf.obj                     = zeros( Float64, nrow(rdf) )                         #objective value
    rdf.opt                     = zeros( Int, nrow(rdf) )                             #optimal ? false, true
    rdf.ones                    = ones( Int, nrow(rdf) )                              #dummy ones
    rdf.SFv                     = Vector{Vector{Int}}(undef, nrow(rdf))               #Followers Seed Set (converted from String to Vector)
    rdf.SLv                     = Vector{Vector{Int}}(undef, nrow(rdf))               #Followers Seed Set (converted from String to Vector)
    rdf.rowID                   = zeros( Int, nrow(rdf) )                             #dummy row ID
    rdf.mtypeString             = Vector{String}(undef, nrow(rdf))
    rdf.obj_resistant_F25       = zeros( Float64, nrow(rdf) )                         #objective value for Fractional model
    rdf.obj_resistant_L25       = zeros( Float64, nrow(rdf) )                         #objective value for Fractional model
    rdf.obj_resistant_F50       = zeros( Float64, nrow(rdf) )                         #objective value for Fractional model
    rdf.obj_resistant_L50       = zeros( Float64, nrow(rdf) )                         #objective value for Fractional model
    rdf.obj_resistant_F75       = zeros( Float64, nrow(rdf) )                         #objective value for Fractional model
    rdf.obj_resistant_L75       = zeros( Float64, nrow(rdf) )                         #objective value for Fractional model
    rdf.obj_total_F             = zeros( Float64, nrow(rdf) )                         #objective value for total views
    rdf.obj_total_L             = zeros( Float64, nrow(rdf) )                         #objective value for total views
    rdf.obj_organic_F           = zeros( Float64, nrow(rdf) )                         #objective value for total views
    rdf.obj_organic_L           = zeros( Float64, nrow(rdf) )                         #objective value for total views
    rdf.obj_forwards_F          = zeros( Float64, nrow(rdf) )                        #objective value for Kempe
    rdf.obj_forwards_L          = zeros( Float64, nrow(rdf) )                         #objective value for Kempe
    rdf.setSimOverModels        = zeros( Float64, nrow(rdf) )                         #seed set similarity over different models
    rdf.wiener_index_f          = zeros( Float64, nrow(rdf) )                         #Wiener Index forwarding
    rdf.wiener_index_v          = zeros( Float64, nrow(rdf) )                         #Wiener index viewing
    rdf.nTotalViewsF            = Vector{Vector{Float64}}(undef, nrow(rdf))
    rdf.nOrganicViewsF          = Vector{Vector{Float64}}(undef, nrow(rdf))
    rdf.nForwardsF              = Vector{Vector{Float64}}(undef, nrow(rdf))
    rdf.nTotalViewsL            = Vector{Vector{Float64}}(undef, nrow(rdf))
    rdf.nOrganicViewsL          = Vector{Vector{Float64}}(undef, nrow(rdf))
    rdf.nForwardsL              = Vector{Vector{Float64}}(undef, nrow(rdf))
    for l=1:nrow(rdf)
        rdf.rowID[l] = l
        rdf.nTotalViewsF[l]     = Float64[]
        rdf.nOrganicViewsF[l]   = Float64[]
        rdf.nForwardsF[l]       = Float64[]
        rdf.nTotalViewsL[l]     = Float64[]
        rdf.nOrganicViewsL[l]   = Float64[]
        rdf.nForwardsL[l]       = Float64[]
        rdf.SFv[l]              = Int[]
        rdf.SLv[l]              = Int[]
        rdf.mtype[l] == 1 ? rdf.mtypeString[l] = "F"     : nothing
        rdf.mtype[l] == 2 ? rdf.mtypeString[l] = "T"     : nothing
        rdf.mtype[l] == 3 ? rdf.mtypeString[l] = "O"     : nothing
        rdf.mtype[l] == 4 ? rdf.mtypeString[l] = "OA"     : nothing
        if rdf.mtype[l] == 5
            rdf.fractionOfResistantNodes[l] == 0.25 ? rdf.mtypeString[l] = "R25" : nothing
            rdf.fractionOfResistantNodes[l] == 0.5  ? rdf.mtypeString[l] = "R50" : nothing
            rdf.fractionOfResistantNodes[l] == 0.75 ? rdf.mtypeString[l] = "R75" : nothing
            if rdf.fractionOfResistantNodes[l] ∉ [0.25, 0.5, 0.75]
                rdf.mtypeString[l] = "R"
            end
        end
    end

    return nothing
end

function convertLeaderUtilitiesToFactors( rdf::DataFrame )
    rdf.leaderUtilityFactor = zeros( Int, nrow(rdf) )
    for l in 1:nrow(rdf)
        rdf.leaderUtility[l] == 0   ? rdf.leaderUtilityFactor[l] = 0 : nothing
        rdf.leaderUtility[l] == 1   ? rdf.leaderUtilityFactor[l] = 1 : nothing
        rdf.leaderUtility[l] == 8   ? rdf.leaderUtilityFactor[l] = 2 : nothing
        rdf.leaderUtility[l] == 20  ? rdf.leaderUtilityFactor[l] = 3 : nothing
        rdf.leaderUtility[l] == 55  ? rdf.leaderUtilityFactor[l] = 4 : nothing
        rdf.leaderUtility[l] == 150 ? rdf.leaderUtilityFactor[l] = 5 : nothing
    end
    return
end


#NOTE: fits for CIMP
function addEmptyDataColumnsAfterEvaluationOfObjectiveFunctions!( rdf::DataFrame )
    rdf.loss_resistant_F25       = zeros( Float64, nrow(rdf) )                         #objective value for Fractional model
    rdf.loss_resistant_L25       = zeros( Float64, nrow(rdf) )                         #objective value for Fractional model
    rdf.loss_resistant_F50       = zeros( Float64, nrow(rdf) )                         #objective value for Fractional model
    rdf.loss_resistant_L50       = zeros( Float64, nrow(rdf) )                         #objective value for Fractional model
    rdf.loss_resistant_F75       = zeros( Float64, nrow(rdf) )                         #objective value for Fractional model
    rdf.loss_resistant_L75       = zeros( Float64, nrow(rdf) )                         #objective value for Fractional model
    rdf.loss_total_F             = zeros( Float64, nrow(rdf) )                         #objective value for total views
    rdf.loss_total_L             = zeros( Float64, nrow(rdf) )                         #objective value for total views
    rdf.loss_organic_F           = zeros( Float64, nrow(rdf) )                         #objective value for total views
    rdf.loss_organic_L           = zeros( Float64, nrow(rdf) )                         #objective value for total views
    rdf.loss_forwards_F          = zeros( Float64, nrow(rdf) )                         #objective value for total views
    rdf.loss_forwards_L           = zeros( Float64, nrow(rdf) )                         #objective value for total views
    rdf.nTotalViewsFv            = Vector{Vector{Float64}}(undef, nrow(rdf))
    rdf.nOrganicViewsFv          = Vector{Vector{Float64}}(undef, nrow(rdf))
    rdf.nForwardsFv              = Vector{Vector{Float64}}(undef, nrow(rdf))
    rdf.nTotalViewsLv            = Vector{Vector{Float64}}(undef, nrow(rdf))
    rdf.nOrganicViewsLv          = Vector{Vector{Float64}}(undef, nrow(rdf))
    rdf.nForwardsLv              = Vector{Vector{Float64}}(undef, nrow(rdf))
    rdf.nTotalViewsFmean         = zeros( Float64, nrow(rdf) )
    rdf.nOrganicViewsFmean       = zeros( Float64, nrow(rdf) )
    rdf.nForwardsFmean           = zeros( Float64, nrow(rdf) )
    rdf.nTotalViewsLmean         = zeros( Float64, nrow(rdf) )
    rdf.nOrganicViewsLmean       = zeros( Float64, nrow(rdf) )
    rdf.nForwardsLmean           = zeros( Float64, nrow(rdf) )

    for l=1:nrow(rdf)
        rdf.nTotalViewsFv[l]     = Float64[]
        rdf.nOrganicViewsFv[l]   = Float64[]
        rdf.nForwardsFv[l]       = Float64[]
        rdf.nTotalViewsLv[l]     = Float64[]
        rdf.nOrganicViewsLv[l]   = Float64[]
        rdf.nForwardsLv[l]       = Float64[]
    end
    return
end


#NOTE: fits for CIMP
function renameAndRoundDataFrameWithEvaluatedValues( rdf )
    for l=1:nrow(rdf)
        rdf.mtype[l] == 1 ? rdf.mtypeString[l] = "F"     : nothing
        rdf.mtype[l] == 2 ? rdf.mtypeString[l] = "T"     : nothing
        rdf.mtype[l] == 3 ? rdf.mtypeString[l] = "O"     : nothing
        rdf.mtype[l] == 4 ? rdf.mtypeString[l] = "OA"     : nothing
        if rdf.mtype[l] == 5
            rdf.fractionOfResistantNodes[l] == 0.25 ? rdf.mtypeString[l] = "R25" : nothing
            rdf.fractionOfResistantNodes[l] == 0.5  ? rdf.mtypeString[l] = "R50" : nothing
            rdf.fractionOfResistantNodes[l] == 0.75 ? rdf.mtypeString[l] = "R75" : nothing
        end
    end
    rdf.obj                  .= round.( rdf.obj ; digits=1 )
    rdf.obj_F_eval           .= round.( rdf.obj_F_eval ; digits=1 )
    rdf.obj_total_F          .= round.( rdf.obj_total_F ; digits=1 )
    rdf.obj_organic_F        .= round.( rdf.obj_organic_F ; digits=1 )
    rdf.obj_forwards_F       .= round.( rdf.obj_forwards_F ; digits=1 )
    rdf.obj_resistant_F25    .= round.( rdf.obj_resistant_F25 ; digits=1 )
    rdf.obj_resistant_F50    .= round.( rdf.obj_resistant_F50 ; digits=1 )
    rdf.obj_resistant_F75    .= round.( rdf.obj_resistant_F75 ; digits=1 )
    rdf.obj_L_eval           .= round.( rdf.obj_L_eval ; digits=1 )
    rdf.obj_total_L          .= round.( rdf.obj_total_L ; digits=1 )
    rdf.obj_organic_L        .= round.( rdf.obj_organic_L ; digits=1 )
    rdf.obj_forwards_L       .= round.( rdf.obj_forwards_L ; digits=1 )
    rdf.obj_resistant_L25    .= round.( rdf.obj_resistant_L25 ; digits=1 )
    rdf.obj_resistant_L50    .= round.( rdf.obj_resistant_L50 ; digits=1 )
    rdf.obj_resistant_L75    .= round.( rdf.obj_resistant_L75 ; digits=1 )

    return
end


function roundValuesForTable!( df, values_to_round, ndigits )
    for mycol in values_to_round
        df[!, mycol] .= round.( df[!,mycol] ; digits = ndigits)
    end
    return
end


# COMPUTE / VALUES
function computeValues!( rdf::DataFrame , p::parameters )
    print("compute values...")
    for l=1:size(rdf,1)
        #rdf[!,:colname][l]
        if( rdf[!,:exitflag][l] ∈ [1,2,3] )   #usable vals (optimal, timelim, memlim)
#            rdf[!,:obj][l] = rdf[!,:obj_F_solcheck][l]
            rdf[!,:obj][l] = rdf[!,:obj_F_eval][l]
            rdf[!,:mtype][l] != 3 ? rdf[!,:rt][l] = min( (rdf[!,:rtCPLEX][l] + rdf[!,:rtTMIN][l] + rdf[!,:rtR][l] +rdf[!,:rtBM][l]), p.tlim ) : nothing
            rdf[!,:gap][l] = round( (rdf[!,:bestbound][l] - rdf[!,:obj][l]) / rdf[!,:bestbound][l] ; digits=4 ) * 100
            rdf[!,:exitflag][l] == 1 ? rdf[!,:gap][l] = 0 : nothing
            rdf[!,:exitflag][l] == 1 ? rdf[!,:opt][l] = 1 : nothing
            rdf[!,:gap][l] == 0 ? rdf[!,:opt][l] = 1 : nothing

        else    #memory exceptions etc
            rdf[!,:gap][l] = 100
            rdf[!,:mtype][l] >= 1 ? rdf[!,:rt][l]  = p.tlim : nothing
            rdf[!,:mtype][l] >= 1 ? rdf[!,:obj][l] = 0      : nothing
        end
    end
    rdf[(rdf.aggoutput .== 1), :obj] .= rdf[(rdf.aggoutput .== 1), :obj_F_eval]

    print("done!\n")
    return nothing
end


#parse seedsets to from string to vector
function parseSeedSets!( rdf )
    print("parse seed sets...")
    #convert seed sets from string to vector
    #function to convert strings to be parse-able
    stripChar  = (s, r) -> replace(s, Regex("[$r]") => "")
    clearInt64 = (s)    -> replace(s, "Int64"       => "")

    #convert df columns
    myStringSFs     = clearInt64.(stripChar.(rdf.SF, "][,"))                     # convert string to e.g. "2 4 4 76", and empty Int64[] -> ""
    myStringSLs     = clearInt64.(stripChar.(rdf.SL, "][,"))

    #parse them
    for l=1:size(rdf,1)
        isempty( myStringSFs[l] )     ? rdf.SFv[l]     = Int[] : rdf.SFv[l]     = parse.( Int, split(myStringSFs[l]) )
        isempty( myStringSLs[l] )     ? rdf.SLv[l]     = Int[] : rdf.SLv[l]     = parse.( Int, split(myStringSLs[l]) )
    end

    print("done!\n")
    return nothing
end

#parse seedsets to from string to vector
function parseSeedSetsAgain!( rdf )
    print("parse seed sets again...")
    rdf.SFvv                     = Vector{Vector{Int}}(undef, nrow(rdf))               #Followers Seed Set (converted from String to Vector)
    for l=1:nrow(rdf)
        rdf.SFvv[l]              = Int[]
    end

    #convert seed sets from string to vector
    #function to convert strings to be parse-able
    stripChar  = (s, r) -> replace(s, Regex("[$r]") => "")
    clearInt64 = (s)    -> replace(s, "Int64"       => "")

    #convert df columns
    myStringSFs     = clearInt64.(stripChar.(rdf.SF, "][,"))                     # convert string to e.g. "2 4 4 76", and empty Int64[] -> ""
#    myStringSLs     = clearInt64.(stripChar.(rdf.SL, "][,"))

    #parse them
    for l=1:size(rdf,1)
        isempty( myStringSFs[l] )     ? rdf.SFvv[l]     = Int[] : rdf.SFvv[l]     = parse.( Int, split(myStringSFs[l]) )
#        isempty( myStringSLs[l] )     ? rdf.SLvv[l]     = Int[] : rdf.SLvv[l]     = parse.( Int, split(myStringSLs[l]) )
    end

    print("done!\n")
    return nothing
end


#parse seedsets to from string to vector #NOTE: fits for CIMP
function parseNumberOfViewings!( rdf )
    print("parse seed sets viewings...")
    #convert seed sets from string to vector
    #function to convert strings to be parse-able
    stripChar  = (s, r) -> replace(s, Regex("[$r]") => "")
    clearFloat64 = (s)    -> replace(s, "Float64"       => "")

    #convert df columns
    myStringTotalViewsF     = clearFloat64.(stripChar.(rdf.nTotalViewsF, "][,"))                     # convert string to e.g. "2 4 4 76", and empty Int64[] -> ""
    myStringOrganicViewsF   = clearFloat64.(stripChar.(rdf.nOrganicViewsF, "][,"))
    myStringForwardsF       = clearFloat64.(stripChar.(rdf.nForwardsF, "][,"))
    myStringTotalViewsL     = clearFloat64.(stripChar.(rdf.nTotalViewsL, "][,"))                     # convert string to e.g. "2 4 4 76", and empty Int64[] -> ""
    myStringOrganicViewsL   = clearFloat64.(stripChar.(rdf.nOrganicViewsL, "][,"))
    myStringForwardsL       = clearFloat64.(stripChar.(rdf.nForwardsL, "][,"))

    #parse them
    for l=1:size(rdf,1)
        isempty( myStringTotalViewsF[l] )      ? rdf.nTotalViewsFv[l]       = Float64[] : rdf.nTotalViewsFv[l]     = parse.( Float64, split(myStringTotalViewsF[l]) )
        isempty( myStringOrganicViewsF[l] )    ? rdf.nOrganicViewsFv[l]     = Float64[] : rdf.nOrganicViewsFv[l]   = parse.( Float64, split(myStringOrganicViewsF[l]) )
        isempty( myStringForwardsF[l] )        ? rdf.nForwardsFv[l]         = Float64[] : rdf.nForwardsFv[l]       = parse.( Float64, split(myStringForwardsF[l]) )
        isempty( myStringTotalViewsL[l] )      ? rdf.nTotalViewsLv[l]       = Float64[] : rdf.nTotalViewsLv[l]     = parse.( Float64, split(myStringTotalViewsL[l]) )
        isempty( myStringOrganicViewsL[l] )    ? rdf.nOrganicViewsLv[l]     = Float64[] : rdf.nOrganicViewsLv[l]   = parse.( Float64, split(myStringOrganicViewsL[l]) )
        isempty( myStringForwardsL[l] )        ? rdf.nForwardsLv[l]         = Float64[] : rdf.nForwardsLv[l]       = parse.( Float64, split(myStringForwardsL[l]) )
    end

    print("done!\n")
    return nothing
end


#my col is the Ranking tuple of interest
function parseSeedSetRanks!( rdf, mycol::Symbol )
    print("parse seed set ranks...")
    stripChar       = (s, r) -> replace(s, Regex("[$r]") => "")
    clearTupleName  = (s)    -> replace(s, "Tuple{Int64Int64}"  => "")

    #add new column

    #convert data to string
    myStringTuples     = clearTupleName.(stripChar.(rdf[!,mycol], "][(),"))

    #parse data
    parsedTupleVector = Vector{Vector{Tuple{Int,Int}}}(undef, nrow(rdf))
    for l=1:size(rdf,1)
        parsedTupleVector[l] = Tuple{Int,Int}[]
        if rdf.aggoutput[l] == 1
            my_vals = parse.( Int, split(myStringTuples[l]) )
#            println(my_vals)
            for i in 1:2:length(my_vals)
                node = my_vals[i]
                rank = my_vals[i+1]
                push!(parsedTupleVector[l], (node, rank))
            end
        end
    end

    #add column
    #insert column of interest
    tempColumnName = string(mycol) * "Vec"
    insertcols!( rdf, ncol(rdf)+1, Symbol(tempColumnName) => parsedTupleVector  )

    print("done!\n")
    return nothing
end



function computeRelativeLossesWRTLeaderUtilityAndSetSimilarity( rdf::DataFrame )
    println("computeRelativeLosses w.r.t leader utility....")
    rdf.loss_leaderUtility       = zeros( Float64, nrow(rdf) )
    rdf.setSim_leaderUtility     = zeros( Float64, nrow(rdf) )
    rdf.avg_nodeFrequency        = zeros( Float64, nrow(rdf) )      #counts how often a seed node appears over all utilities

    for subdf_file in groupby(rdf, :ifile)
        for subdf_nsc in groupby(subdf_file, :nsc)
            for subdf_kF in groupby( subdf_nsc, :kF)
                for subdf_kL in groupby(subdf_kF, :kL)
                    for subdf_mtype in groupby(subdf_kL, :mtype)
                        for subdf_fr in groupby(subdf_mtype, :fractionOfResistantNodes)
                            mydf = subdf_fr
                            df_lu = mydf |> @filter(_.aggoutput == 1 ) |> DataFrame
                            base_value = df_lu.obj_F_eval[1]
                            myIntersection = df_lu[!,:SFv][1]
                            mySeedSetCollection = Int[]
                            for l in 1:nrow(df_lu)
                                loss    = (df_lu.obj_F_eval[l] / base_value) * 100
                                rowID   = df_lu.rowID[l]
                                rdf.loss_leaderUtility[rowID] = loss
                                myIntersection = intersect( myIntersection, df_lu[!,:SFv][l] )
                                append!(mySeedSetCollection, df_lu[!,:SFv][l])
                            end
                            intersection_ratio = length( myIntersection ) * 100 / df_lu.kF[1]
                            myFrequencyDict    = collect(StatsBase.countmap(mySeedSetCollection))
                            avg_nodeFrequency  = 0.0
                            for i in 1:length(myFrequencyDict)
                                avg_nodeFrequency += myFrequencyDict[i][2]
                            end
                            avg_nodeFrequency = ((avg_nodeFrequency / length(myFrequencyDict)) / 6) * 100
                            for l in 1:nrow(df_lu)
                                rowID   = df_lu.rowID[l]
                                rdf.setSim_leaderUtility[rowID] = intersection_ratio
                                rdf.avg_nodeFrequency[rowID]    = avg_nodeFrequency
                            end
                        end
                    end
                end
            end
        end
    end

    println("done!")
    return
end


#compute seed set similarity over all SAA iterations
function computeSetSimilarityOverSAAIterations!( rdf )
    print("compute set similarity...")
    rdf.setSim  = zeros( Float64, nrow(rdf) )                       #set similarity over all SSAItr (the aggregared file will contain the result)

    #compute set similarity
    for subdf_file in groupby(rdf, :ifile)
        for subdf_nsc in groupby(subdf_file, :nsc)
            for subdf_kF in groupby( subdf_nsc, :kF)
                for subdf_kL in groupby(subdf_kF, :kL)
                    for subdf_mtype in groupby(subdf_kL, :mtype)
                        for subdf_fr in groupby(subdf_mtype, :fractionOfResistantNodes)
                            mydf = subdf_fr
                            if nrow(mydf) == 11       #number of SAA iterations + aggregated file
                                df_SAA = mydf |> @filter(_.aggoutput == 0 ) |> DataFrame
                                myIntersection = df_SAA[!,:SFv][1]
                                for l in 1:nrow(df_SAA)
                                    myIntersection = intersect( myIntersection, df_SAA[!,:SFv][l] )
                                end
                                intersection_ratio = length( myIntersection ) * 100 / df_SAA.kF[1]      #in percent

                                df_AGG = mydf |> @filter(_.aggoutput == 1 ) |> DataFrame
                                rowID_in_rdf = df_AGG.rowID[1]
                                rdf.setSim[rowID_in_rdf] = intersection_ratio
                            else
                                my_nrow = nrow(mydf)
                                @warn "function computeSetSimilarityOverSAAIterations: more than 11 files to process ($my_nrow)"
                            end
                        end
                    end
                end
            end
        end
    end
    print("done!\n")

    return
end


function computeSetSimilarityRmodels!( df )
#    df.setSimR  = zeros( Float64, nrow(df) )

    for subdf_ifile in groupby(df, :ifile)
        println("")
        for subdf_kF in groupby(subdf_ifile, :kF)
            mydf = subdf_kF
            myIntersection = mydf[!,:SFvv][1]
            myUnion = Int[]
            for l in 1:nrow(mydf)
                myIntersection = intersect( myIntersection, mydf[!,:SFvv][l] )
                myUnion = union(myUnion, mydf[!,:SFvv][l])
            end
            intersection_ratio = length( myIntersection ) * 100 / length(myUnion)
            println("kF: ", mydf.kF[1], "\tSim: ", intersection_ratio)

            for l in 1:nrow(mydf)
                mydf.setSimR[l] = intersection_ratio
            end
        end
    end

    return
end


function computeSampleStability!( rdf )
    print("compute sample stability...")
    rdf.inStability    = zeros( Float64, nrow(rdf) )                            #in-sample stability
    rdf.outStability   = zeros( Float64, nrow(rdf) )                            #out of sample stability

    #compute set similarity
    for subdf_file in groupby(rdf, :ifile)
        for subdf_nsc in groupby(subdf_file, :nsc)
            for subdf_kF in groupby( subdf_nsc, :kF)
                for subdf_kL in groupby(subdf_kF, :kL)
                    for subdf_mtype in groupby(subdf_kL, :mtype)
                        for subdf_fr in groupby(subdf_mtype, :fractionOfResistantNodes)
                            mydf = subdf_fr
                            if nrow(mydf) == 11       #number of SAA iterations + aggregated file
                                df_AGG = mydf |> @filter(_.aggoutput == 1 ) |> DataFrame
                                rowID_in_rdf = df_AGG.rowID[1]

                                df_stab = mydf |> @filter( (_.aggoutput == 0) && (_.exitflag == 1) ) |> DataFrame
                                my_instability_vals  = Float64[]
                                my_outstability_vals = Float64[]
                                if nrow(df_stab) >= 2
                                    for k in 1:nrow(df_stab)
                                        push!(my_instability_vals,  df_stab.obj_F_solcheck[k])
                                        push!(my_outstability_vals, df_stab.obj_F_eval[k])
                                    end
                                    sort!(my_instability_vals)
                                    sort!(my_outstability_vals)

                                    counter = 1
                                    avg_instability = 0.0
                                    avg_outstability = 0.0
                                    while( length(my_instability_vals) != 1)
                                        lowest_val_in  = popfirst!(my_instability_vals)
                                        lowest_val_out = popfirst!(my_outstability_vals)
                                        for k in 1:length(my_instability_vals)
                                            avg_instability  += (1 - lowest_val_in / my_instability_vals[k]) * 100
                                            avg_outstability += (1 - lowest_val_out / my_outstability_vals[k]) * 100
                                            counter += 1
                                        end
                                    end
                                    avg_instability  = avg_instability  / counter
                                    avg_outstability = avg_outstability / counter
                                    rdf.inStability[rowID_in_rdf] = avg_instability
                                    rdf.outStability[rowID_in_rdf] = avg_outstability
                                end
                            else
                                #println(mydf.ifile[1], "\tw=", mydf.nsc[1], "\tkL=", mydf.kL[1], "\tkF=", mydf.kF[1], "\tNOT COMPLETE!! (n=", nrow(mydf), ")")
                            end
                        end
                    end
                end
            end
        end
    end
    print("done!\n")

    return
end


function compareObjectiveValuesRelativeToScenarioSize( rdf )
    println("compare objective values...")

    #compute set similarity
    df_AGG = rdf |> @filter(_.aggoutput == 1 ) |> DataFrame
    for subdf_file in groupby(df_AGG, :ifile)
        print(subdf_file.ifile[1], "\n")
        for subdf_kF in groupby( subdf_file, :kF)
            print("\tkF: ", subdf_kF.kF[1])
            for subdf_nsc in groupby(subdf_kF, :nsc)
                print("\t", subdf_nsc.obj[1])
            end
            print("\n")
        end
    end
    print("done!\n")

    return
end


function getStackedDataFrameForComparisonOfObjectiveValuesWithHeursitic( rdf )
    df = rdf |> @filter(_.aggoutput == 1 ) |> DataFrame
    my_heuristics = [:obj_F_eval, :obj_FMAR_eval, :obj_Fdc_eval, :obj_Fed_eval, :obj_Fbc_eval, :obj_Fpr_eval, :obj_Ftr_eval, :obj_Frm_eval]
    my_heuristic_names = ["BEN", "MAR", "DC", "EC", "BC", "PR", "TR", "RM"]

    #add column used as category for boxplots
    tempStringVector  = String[]
    for i in 1:nrow(df)
        push!(tempStringVector, "BEN")
    end
    df.methodname = tempStringVector
    df.obj_ratio    = zeros( Float64, nrow(rdf) )

    stacked_df = stack(df, my_heuristics)
    for l in 1:nrow(stacked_df)
        idx_of_heuristic = findfirst((x -> x == stacked_df.variable[l]), my_heuristics)           #index of variable in my_heuristic Vector
        stacked_df[!,:methodname][l]   = my_heuristic_names[idx_of_heuristic]
        stacked_df[!,:obj_ratio][l]    = (stacked_df.value[l] / stacked_df[!,:obj][l]) * 100
        stacked_df[!,:obj][l]          = stacked_df.value[l]
    end

    stacked_df = stacked_df |> @filter(_.methodname != "BEN" ) |> DataFrame
    filter!(row -> !(row.methodname == "RM" && row.ifile == "msg-college"),  stacked_df)
    filter!(row -> !(row.methodname == "RM" && row.ifile == "msg-email-eu"),  stacked_df)
    filter!(row -> !(row.methodname == "RM" && row.ifile == "soc-advogato"),  stacked_df)
    filter!(row -> !(row.methodname == "RM" && row.ifile == "soc-anybeat"),  stacked_df)
    filter!(row -> !(row.methodname == "TR" && row.ifile == "msg-college"),  stacked_df)
    filter!(row -> !(row.methodname == "TR" && row.ifile == "msg-email-eu"),  stacked_df)
    filter!(row -> !(row.methodname == "TR" && row.ifile == "soc-advogato"),  stacked_df)
    filter!(row -> !(row.methodname == "TR" && row.ifile == "soc-anybeat"),  stacked_df)

    return stacked_df
end


#NOTE: does not respect leader views
function getStackedDataFrameForComparisonOfDifferentMethods( rdf )
    df = rdf |> @filter(_.aggoutput == 1 ) |> DataFrame
    my_methods = [:loss_forwards_F, :loss_organic_F, :loss_total_F, :loss_resistant_F25, :loss_resistant_F50, :loss_resistant_F75]
    my_method_names = ["F", "O", "T", "R25", "R50", "R75"]

    #add column used as category for boxplots
    tempStringVector  = String[]
    for i in 1:nrow(df)
        push!(tempStringVector, "F")
    end
    df.methodname = tempStringVector
    df.losses    = zeros( Float64, nrow(df) )

    stacked_df = stack(df, my_methods)
#    print(stacked_df.variable)
    for l in 1:nrow(stacked_df)
        idx_of_methods = findfirst((x -> x == stacked_df.variable[l]), my_methods)           #index of variable in my_heuristic Vector
#        println(idx_of_methods)
        stacked_df[!,:methodname][l]   = my_method_names[idx_of_methods]
        stacked_df[!,:losses][l]       = stacked_df.value[l]                                #value is the stacked value
    end

    return stacked_df
end


function computeSeedSetSimilarityOverDifferentModels!( rdf )
    print("compute seed set similarity over different models...")
    df_AGG = rdf |> @filter(_.aggoutput == 1 ) |> DataFrame

    if length(unique(df_AGG.nsc)) > 1
        @warn "function computeSeedSetSimilarityOverDifferentModels: more than one type of |Ω| to consider"
    end

    for subdf_file in groupby(df_AGG, :ifile)
        for subdf_kL in groupby(subdf_file, :kL)
            for subdf_kF in groupby(subdf_kL, :kF)
                #mydf = subdf_kF
                mydf = subdf_kF
                kL = mydf.kL[1]
                kF = mydf.kF[1]
                intersection = mydf.SFv[1]
                for i in 2:size(mydf,1)
                    intersection = intersect(intersection, mydf.SFv[i])
                end
                setSimilarity = length(intersection) / kF
            #    println(setSimilarity)
                for i in 1:size(mydf,1)
                    rowID = mydf.rowID[i]
                    rdf.setSimOverModels[rowID] = setSimilarity
                end
            end
        end
    end

    print("done!\n")
    return
end

# reduce objective value of classical approach by |F| (for fair)
function reduceForwardingObjectiveValueBySeedSetSize!( rdf::DataFrame )
    for l in 1:nrow(rdf)
        if rdf.obj_forwards_F[l] >= 0
            rdf.obj_forwards_F[l] -= rdf.kF[l]
            if rdf.mtypeString[l] == "F"
                rdf.obj[l] -= rdf.kF[l]
                rdf.obj_F_eval[l] -= rdf.kF[l]
            end
        end
    end
    return
end


#NOTE: fits for CIMP
function computeRelativeLossesWithRespectToDifferentModels!( rdf )
    df_AGG = rdf |> @filter(_.aggoutput == 1 ) |> DataFrame

    for subdf_file in groupby(df_AGG, :ifile)
        #group by kL
        for subdf_kL in groupby(subdf_file, :kL)
            #groupby kF
            for subdf_kF in groupby(subdf_kL, :kF)
                mydf = subdf_kF
#                sort!(mydf, [:mtypeString, :kF])
#                mydf = mydf[[:ifile, :kF, :kL, :mtypeString, :obj]]
                idx_forwards    = findall(x->x=="F", mydf.mtypeString)[1]
                idx_organic     = findall(x->x=="O", mydf.mtypeString)[1]
                idx_total       = findall(x->x=="T", mydf.mtypeString)[1]
                idx_R25         = findall(x->x=="R25", mydf.mtypeString)[1]
                idx_R50         = findall(x->x=="R50", mydf.mtypeString)[1]
                idx_R75         = findall(x->x=="R75", mydf.mtypeString)[1]
                for i in 1:nrow(mydf)
                    rowID       = mydf.rowID[i]
                    obj_val     = mydf.obj[i]
                    obj_val_L   = mydf.obj_L_eval[i]
                    rdf.loss_total_F[rowID]             = (rdf.obj_total_F[rowID] / mydf.obj[idx_total]) * 100
                    rdf.loss_organic_F[rowID]           = (rdf.obj_organic_F[rowID] / mydf.obj[idx_organic]) * 100
                    rdf.loss_forwards_F[rowID]          = (rdf.obj_forwards_F[rowID] / mydf.obj[idx_forwards]) * 100
                    rdf.loss_resistant_F25[rowID]       = (rdf.obj_resistant_F25[rowID] / mydf.obj[idx_R25]) * 100
                    rdf.loss_resistant_F50[rowID]       = (rdf.obj_resistant_F50[rowID] / mydf.obj[idx_R50]) * 100
                    rdf.loss_resistant_F75[rowID]       = (rdf.obj_resistant_F75[rowID] / mydf.obj[idx_R75]) * 100
                    if mydf.kL[i] != 0
                        rdf.loss_total_L[rowID]             = (rdf.obj_total_L[rowID] / mydf.obj_L_eval[idx_total]) * 100
                        rdf.loss_organic_L[rowID]           = (rdf.obj_organic_L[rowID] / mydf.obj_L_eval[idx_organic]) * 100
                        rdf.loss_forwards_L[rowID]          = (rdf.obj_forwards_L[rowID] / mydf.obj_L_eval[idx_forwards]) * 100
                        rdf.loss_resistant_L25[rowID]       = (rdf.obj_resistant_L25[rowID] / mydf.obj_L_eval[idx_R25]) * 100
                        rdf.loss_resistant_L50[rowID]       = (rdf.obj_resistant_L50[rowID] / mydf.obj_L_eval[idx_R50]) * 100
                        rdf.loss_resistant_L75[rowID]       = (rdf.obj_resistant_L75[rowID] / mydf.obj_L_eval[idx_R75]) * 100
                    end
                end
            end
        end
    end

    rdf.loss_total_F          .= round.( rdf.loss_total_F ; digits=1 )
    rdf.loss_organic_F        .= round.( rdf.loss_organic_F ; digits=1 )
    rdf.loss_forwards_F       .= round.( rdf.loss_forwards_F ; digits=1 )
    rdf.loss_resistant_F25    .= round.( rdf.loss_resistant_F25 ; digits=1 )
    rdf.loss_resistant_F50    .= round.( rdf.loss_resistant_F50 ; digits=1 )
    rdf.loss_resistant_F75    .= round.( rdf.loss_resistant_F75 ; digits=1 )
    rdf.loss_total_L          .= round.( rdf.loss_total_L ; digits=1 )
    rdf.loss_organic_L        .= round.( rdf.loss_organic_L ; digits=1 )
    rdf.loss_forwards_L       .= round.( rdf.loss_forwards_L ; digits=1 )
    rdf.loss_resistant_L25    .= round.( rdf.loss_resistant_L25 ; digits=1 )
    rdf.loss_resistant_L50    .= round.( rdf.loss_resistant_L50 ; digits=1 )
    rdf.loss_resistant_L75    .= round.( rdf.loss_resistant_L75 ; digits=1 )

    return
end

#NOTE: fits for CIMP
function cleanOutlier!(rdf)
    #println("in function")
    for l in 1:nrow(rdf)
        rdf.loss_total_F[l]          = min(rdf.loss_total_F[l], 100.0)
        rdf.loss_organic_F[l]        = min(rdf.loss_organic_F[l], 100.0)
        rdf.loss_forwards_F[l]       = min(rdf.loss_forwards_F[l], 100.0)
        rdf.loss_resistant_F25[l]    = min(rdf.loss_resistant_F25[l] , 100.0)
        rdf.loss_resistant_F50[l]    = min(rdf.loss_resistant_F50[l], 100.0)
        rdf.loss_resistant_F75[l]    = min(rdf.loss_resistant_F75[l], 100.0)
        rdf.loss_total_L[l]          = min(rdf.loss_total_L[l], 100.0)
        rdf.loss_organic_L[l]        = min(rdf.loss_organic_L[l], 100.0)
        rdf.loss_forwards_L[l]       = min(rdf.loss_forwards_L[l], 100.0)
        rdf.loss_resistant_L25[l]    = min(rdf.loss_resistant_L25[l], 100.0)
        rdf.loss_resistant_L50[l]    = min(rdf.loss_resistant_L50[l], 100.0)
        rdf.loss_resistant_L75[l]    = min(rdf.loss_resistant_L75[l], 100.0)
    end
    return
end


function removeNonOptimalSolutions(rdf)
    println("removing non optimal instances...")
    line_to_lookup = 0
    start_of_rows_to_delete = Int[]
    dummy = 0
    for l in 1:nrow(rdf)
        if rdf.aggoutput[l] == 1
            line_to_lookup = l + rdf.best_ss_id[l]
            dummy = l
        end
        if l == line_to_lookup
            #println("looking up line $l")
            if rdf.exitflag[l] != 1
                #@warn "not optimal in line $l"
                println(rdf.ifile[l], "\t", rdf.mtypeString[l], "\tkL=", rdf.kL[l], "\tkF=", rdf.kF[l])
                push!(start_of_rows_to_delete, dummy)
            end
        end
    end
    sort!(start_of_rows_to_delete, rev=true)
#    println(start_of_rows_to_delete)

    clean_df = deepcopy(rdf)
    for k in start_of_rows_to_delete
        for i in k+10:-1:k
            clean_df = clean_df[setdiff(1:end, i), :]
        end
    end

    println("done!")
    return clean_df
end


function getJoinedDataFrameWithRunTimes( df_AGG )
    df_rt = df_AGG[:, [:ifile, :kF, :mtypeString, :rtCPLEX]]
    df_joined = DataFrame()
    counter = 0
    for subdf_method in groupby(df_rt, :mtypeString)
        mydf = DataFrame(subdf_method)                      #converts SubDataFrame to DataFrame
        if counter == 0
            df_joined = mydf
            counter += 1
        else
            df_joined = join(df_joined, mydf, on = [:ifile, :kF], kind = :inner, makeunique=true)
        end
    end

    return df_joined
end


function getJoinedDataFrameWithRunTimesCIMP( df_AGG )
    df_rt = df_AGG[:, [:ifile, :kF, :leaderUtility, :rtCPLEX]]
    df_joined = DataFrame()
    counter = 0
    for subdf_lu in groupby(df_rt, :leaderUtility)
        mydf = DataFrame(subdf_lu)                      #converts SubDataFrame to DataFrame
        if counter == 0
            df_joined = mydf
            counter += 1
        else
            df_joined = join(df_joined, mydf, on = [:ifile, :kF], kind = :inner, makeunique=true)
        end
    end

    df_obj = df_AGG[:, [:ifile, :kF, :leaderUtility, :obj_F_eval]]
    for subdf_lu in groupby(df_obj, :leaderUtility)
        mydf = DataFrame(subdf_lu)                      #converts SubDataFrame to DataFrame
        if counter == 0
            df_joined = mydf
            counter += 1
        else
            df_joined = join(df_joined, mydf, on = [:ifile, :kF], kind = :inner, makeunique=true)
        end
    end

    return df_joined
end


#NOTE: fits for CIMP
function computeMeanNumbersOfViewsAndForwards!(rdf)
    for l in 1:nrow(rdf)
        isempty(rdf.nTotalViewsFv[l])       ? nothing : rdf.nTotalViewsFmean[l]     = mean(rdf.nTotalViewsFv[l])
        isempty(rdf.nOrganicViewsFv[l])     ? nothing : rdf.nOrganicViewsFmean[l]   = mean(rdf.nOrganicViewsFv[l])
        isempty(rdf.nForwardsFv[l])         ? nothing : rdf.nForwardsFmean[l]       = mean(rdf.nForwardsFv[l])
        isempty(rdf.nTotalViewsLv[l])       ? nothing : rdf.nTotalViewsLmean[l]     = mean(rdf.nTotalViewsFv[l])
        isempty(rdf.nOrganicViewsLv[l])     ? nothing : rdf.nOrganicViewsLmean[l]   = mean(rdf.nOrganicViewsLv[l])
        isempty(rdf.nForwardsLv[l])         ? nothing : rdf.nForwardsLmean[l]       = mean(rdf.nForwardsLv[l])
    end

    return
end


function roundValuesOfInterest!(rdf)
    rdf.rt                   .= round.( Int, rdf.rt )
    rdf.rtHMAR               .= round.( Int, rdf.rtHMAR)
    rdf.gap_approx_rel       .= round.( rdf.gap_approx_rel ; digits=2 )
    rdf.setSimOverModels     .= round.( rdf.setSimOverModels ; digits=2 ) .* 100
    rdf.setSim               .= round.( rdf.setSim ; digits=2 ) .* 100

    for l in 1:nrow(rdf)
        rdf.nTotalViewsF[l]     = round.(rdf.nTotalViewsF[l], digits=4)
        rdf.nTotalViewsL[l]     = round.(rdf.nTotalViewsL[l], digits=4)
        rdf.nOrganicViewsF[l]   = round.(rdf.nOrganicViewsF[l], digits=4)
        rdf.nOrganicViewsL[l]   = round.(rdf.nOrganicViewsL[l], digits=4)
        rdf.nForwardsF[l]       = round.(rdf.nForwardsF[l], digits=4)
        rdf.nForwardsL[l]       = round.(rdf.nForwardsL[l], digits=4)
    end

    return
end

