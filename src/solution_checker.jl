# CHECK SOLUTION
function check_solution!( inst::instance, res::results, b::bfsdata, model::Model )
    print("check solution...")
    res.sol_valid  = 1

    # retrieve solution
    if inst.kL == 0  #imp
        x_val = annotate( value.( model[:x] ) )
        @assert( length( x_val ) == inst.nnodes, "|x| != |V|" )

        # get seed set
        for i ∈ 1:inst.nnodes
            if( x_val[i] > 0.99999 )
                push!( inst.SF, i )
            end
        end
    else            #cimp
        y_val = annotate( value.( model[:y] ) )
        @assert( length( y_val ) == inst.nnodes, "|y| != |V|" )

        # get seed set
        for i ∈ 1:inst.nnodes
            if( y_val[i] > 0.99999 )
                push!( inst.SF, i )
            end
        end
    end
    @assert( length(inst.SF) <= inst.kF, "|S^F| > kF" )

    # compute objective values
    println("\nPropagate from F: ", inst.SF)

    if inst.kL == 0         #imp
        res.obj_F_solcheck     = computeObjectiveValuesFromSetIMP( inst, b, inst.SF )
    else                    #cimp
        res.obj_F_solcheck     = computeObjectiveValuesFromSetCIMP( inst, b, inst.SF )
    end
    res.fcut_error     = (1 - res.obj_F / res.obj_F_solcheck) * 100

    # check validity of solution
    if (res.cpxStatus == "OPTIMAL") || (res.cpxStatus == "ALMOST_OPTIMAL")
        if abs( res.fcut_error ) > 1 #percent
            res.sol_valid = 0
        end
    end 

    # print output
    if( res.sol_valid == 1 )
        @info "valid!"
        print("obj val (Sol_checker): ", res.obj_F_solcheck, "\n")
    else
        @warn "NOT valid!!"
        print("obj val (Sol_checker): ", res.obj_F_solcheck, "\n")
    end
    return nothing
end # check_solution


#COMPUTE APPROXIMATION GAP
function computeApproximationGap( inst::instance, res::results, u_vals::Vector{Float64}, bestbound::Vector{Float64} )
    gap::Float64 = 0.0
    α::Float64   = 0.01     #confidence level

    z = quantile( Normal(), 1-α)
    t = quantile( TDist(inst.nSSAItr -1), 1-α )

    #compute lower bound (w.r.t evaluation scenarios)
    avg = mean(u_vals)
    var = 1/(inst.Nsc * (inst.Nsc-1)) * sum( (u_vals[i] - avg)^2 for i in 1:inst.Nsc)
    LN = avg - z * sqrt( var )                      #lower bound on the true objective function

    #compute upper  bound (w.r.t. SAA iterations)
    avg = mean(bestbound)
    var = 1/(inst.nSSAItr * (inst.nSSAItr-1)) * sum( (bestbound[i] - avg)^2 for i in 1:inst.nSSAItr)
    UR = avg + t * sqrt(var)     #upper bounds on the expected value of the best known bounds)

    #compute gap
    gap = ((UR - LN) / UR ) * 100

    println("approx. gap: $gap")

    return gap::Float64
end
