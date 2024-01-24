#****************************************************************
#**** MODELS IMP GENERALIZED BENDERS DECOMPOSITION
#****************************************************************
function getInitialBendersCuts( inst::instance )
    ρ::Matrix{Float64} = zeros( Float64, inst.nnodes, inst.nsc )                #cut coefficients
    α::Vector{Float64} = zeros( Float64, inst.nnodes )                          #the same alphas ∀ ω ∈ Ω

    #compute α[i,w] = α[i] ∀ ω ∈ Ω
    if inst.mtype == 5      #fractional benders
        for i=1:inst.nnodes
            α[i] = 1 / inst.r[i]
        end
    end
    if inst.mtype in [2,3]  #total / organic views
        for i=1:inst.nnodes
            α[i] = 1
        end
    end

    #compute ρ[i,w]
    if inst.mtype in [2,3,5]
        for w=1:inst.nsc
            for i=1:inst.nnodes
                for k ∈ inst.Rout[i,w]
                    for j ∈ inst.lgo[k,w]
                        ρ[i,w] += α[j]
                    end
                end
            end
        end
    end

    if inst.mtype == 1      #forwarders
        for w=1:inst.nsc
            for i=1:inst.nnodes
                for k ∈ inst.Rout[i,w]
                    ρ[i,w] += 1
                end
            end
        end
    end

    return ρ::Matrix{Float64}
end


# MODEL GENERALIZED BENDERS DECOMPOSITION
function buildModel_BEN( inst::instance, params::Dict, model::Model, to::TimerOutput )
    print("build model...")

    # get initial cut coefficients
    ρ::Matrix{Float64} = getInitialBendersCuts( inst )

    #vars
    @variable( model, u[w=1:inst.nsc] >= 0 )                           #scenario contribution
    @variable( model, x[i=1:inst.nnodes], Bin, base_name="x" )         #seed set variables

    #obj
    @objective( model, Max, (1/inst.nsc) * sum( u ) )

    #constraints
    @constraint( model, cardinality, sum( x[i] for i=1:inst.nnodes ) <= inst.kF )            # |F| <= k
    @constraint( model, initCuts[w=1:inst.nsc], u[w] <= sum( ρ[i,w] * x[i] for i=1:inst.nnodes) )   # Benders Cut with x=0 ∀i ∈ V

    print("done! \n")
    return
end



#****************************************************************
#**** GBD CALLBACKS
#****************************************************************
# all models with preprocessing
function computeSubgradient_all_models!( inst::instance, res::results, b::bfsdata, S::Vector{Tnodes}, w::Int, wherefrom::String )

    #compute α_i / φ_i
    #compute α_i / φ_i
    if inst.mtype == 5  #GBD
        for i ∈ 1:inst.nnodes
            if !inst.Singleton[i,w]                     #dual vars do not exist for singletons
                if ( b.xv[i] >= inst.indeg[i,w] )
                    b.phi[i] = inst.r[i] / (inst.indeg[i,w] + inst.r[i])^2
                else
                    b.alpha[i] = inst.r[i] / (b.xv[i] + inst.r[i])^2
                end
            end
        end
        for i ∈ S
            b.phi[i] = 1 / inst.r[i]
            b.alpha[i] = 0.0
        end
    end
    if inst.mtype in [2,3]      #total/organic views
        for i ∈ 1:inst.nnodes
            if !inst.Singleton[i,w]                     #dual vars do not exist for singletons
                if ( b.xv[i] >= inst.indeg[i,w] )
                    b.phi[i] = 1
                else
                    b.alpha[i] = 1
                end
            end
        end
        for i ∈ S
            b.phi[i] = 1
            b.alpha[i] = 0.0
        end
    end

    #compute ρ_i
    if inst.mtype in [2,3,5]
        for i ∈ 1:inst.nnodes
            if !inst.Singleton[i,w]                     #dual vars do not exist for singletons
                for k ∈ inst.Rout[i,w]
                    if !b.xfb[k]                              #β_k = 0 ∀ forwarders
                        for j ∈ inst.lgo[k,w]
                            b.rho[i] += b.alpha[j]
                        end
                    end
                end
                b.rho[i] -= inst.indeg[i,w] * b.phi[i]
            end
        end
    end
    

    if inst.mtype == 1         #forwarding
        #compute β_i / γ_i NOTE: b.alpha is used here instead of b.beta, and γ is never used
        if wherefrom == "LCB"
            for i ∈ 1:inst.nnodes
                if ( b.xfb[i] == 0 )
                    b.alpha[i] = 1
                end
            end
        end
        if wherefrom == "UCB"
            for i ∈ 1:inst.nnodes
                if ( b.xf[i] < 1 )
                    b.alpha[i] = 1
                end
            end
        end

        for i ∈ 1:inst.nnodes
            for j ∈ inst.Rout[i,w]
                b.rho[i] += b.alpha[j]
            end
        end
    end

    return
end


#LAZY BENDERS CALLBACK
function computeLazyObjective!( inst::instance, res::results, b::bfsdata, S::Vector{Tnodes}, w::Int )
    #compute xf, xv
    for i ∈ S
        for k ∈ inst.Rout[i,w]
            if !b.xfb[k]                     #forwarding not marked yet
                b.xfb[k] = true
                for j ∈ inst.lgo[k,w]
                    b.xv[j] += 1.0
                end
            end
        end
    end
    for i ∈ S
        b.xv[i] = 0.0
    end

    #adjust viewings for organic views
    if inst.mtype == 3
        for i in 1:inst.nnodes
            if b.xv[i] > 0
                b.xv[i] = 1.0
            end
        end
    end

    #compute solution of Benders subproblem
    if inst.mtype == 5          #fractional benders
        for i ∈ 1:inst.nnodes
            b.Phi  += b.xv[i] / (b.xv[i] + inst.r[i])
        end
    end
    if inst.mtype in [2,3]      #total/organic views
        for i ∈ 1:inst.nnodes
            b.Phi  += b.xv[i]
        end
    end
    if inst.mtype == 1          #forwarders (Kempe)
        for i ∈ 1:inst.nnodes
            b.Phi  += b.xfb[i]
        end
    end

    return
end


#LAZY CUTS GBD
function startLazyBendersCallback( cb_data::CPLEX.CallbackContext, inst::instance, params::Dict, res::results, b::bfsdata, model::Model )
    relcuttol::Float64 = params["relcuttol"]
    wherefrom::String = "LCB"

    # get current values
    x_val = annotate( callback_value.(Ref(cb_data), model[:x]) )
    u_val = annotate( callback_value.(Ref(cb_data), model[:u]) )

    # set current seed set
    S = Tnodes[]
    sizehint!(S, inst.kF)
    for i ∈ eachindex( x_val )
        if x_val[i] > 0.5
            push!(S, i)
        end
    end

    # add violated cuts
    for w ∈ 1:inst.nsc
        clear!(b)
        computeLazyObjective!( inst, res, b, S, w )

        if (u_val[w] / b.Phi) > relcuttol
            computeSubgradient_all_models!( inst, res, b, S, w, wherefrom )            

            #add cuts
            con = @build_constraint( model[:u][w] <= b.Phi + sum( b.rho[i] * model[:x][i] for i ∈ 1:inst.nnodes ) - sum( b.rho[i] for i ∈ S ) )
            MOI.submit(model, MOI.LazyConstraint(cb_data), con)
            res.nBENCUTS += 1
            res.nLCUTS   += 1
        end
    end # ∀ω ∈ Ω

    return
end

#------- USER CUTS
function computeUserObjective!( inst::instance, res::results, b::bfsdata, x_val::Vector{Float64}, w::Int )
    tempsum::Float64 = 0.0

    #compute xf
    for i ∈ 1:inst.nnodes
        if !inst.Singleton[i,w]                         #all values are zero for singletons
            tempsum = 0.0
            for k ∈ inst.Rin[i,w]
                tempsum += x_val[k]
                tempsum >= 1.0 ? break : nothing
            end
            b.xf[i] = min( 1.0, tempsum )
            if b.xf[i] == 1.0
                b.xfb[i] = true
            end
        end
    end

    #compute xv
    if inst.mtype == 5
        for i ∈ 1:inst.nnodes
            if !inst.Singleton[i,w]                         #all values are zero for singletons
                tempsum = 0.0
                for j ∈ inst.lgi[i,w]
                    tempsum += b.xf[j]
                end
                b.xv[i] = min( inst.indeg[i,w] * (1-x_val[i]), tempsum )
                b.Phi  += b.xv[i] / (b.xv[i] + inst.r[i])
            end
        end
    end
    if inst.mtype == 2      #total views
        for i ∈ 1:inst.nnodes
            if !inst.Singleton[i,w]                         #all values are zero for singletons
                tempsum = 0.0
                for j ∈ inst.lgi[i,w]
                    tempsum += b.xf[j]
                end
                b.xv[i] = min( inst.indeg[i,w] * (1-x_val[i]), tempsum )
                b.Phi  += b.xv[i]
            end
        end
    end
    if inst.mtype == 3      #organic views
        for i ∈ 1:inst.nnodes
            if !inst.Singleton[i,w]                         #all values are zero for singletons
                tempsum = 0.0
                for j ∈ inst.lgi[i,w]
                    tempsum += b.xf[j]
                end
                b.xv[i] = min( inst.indeg[i,w] * (1-x_val[i]), tempsum )
                if b.xv[i] > 1
                    b.xv[i] = 1
                end
                b.Phi  += b.xv[i]
            end
        end
    end
    if inst.mtype == 1      #Kempe
        for i in 1:inst.nnodes
            b.Phi += b.xf[i]
        end
    end

    return
end


#USER BENDERS CALLBACK
function startUserBendersCallback( cb_data::CPLEX.CallbackContext, inst::instance, params::Dict, res::results, b::bfsdata, model::Model, to::TimerOutput )
    eps::Float64 = 0.001
    relcuttol::Float64 = params["relcuttol"]
    wherefrom::String = "UCB"

    x_val = annotate( callback_value.(Ref(cb_data), model[:x]) )
    u_val = annotate( callback_value.(Ref(cb_data), model[:u]) )

    # set current seed set
    S::Vector{Tnodes} = Tnodes[]
    sizehint!(S, inst.kF)
    for i ∈ eachindex( x_val )
        if x_val[i] > 1-eps
            x_val[i] = 1.0
            push!(S, i)
        end
        if x_val[i] < eps
            x_val[i] = 0.0
        end
    end

    # add violated cuts
    for w ∈ 1:inst.nsc
        clear!(b)
        computeUserObjective!( inst, res, b, x_val, w )

        if (u_val[w] / b.Phi) > relcuttol
            computeSubgradient_all_models!( inst, res, b, S, w, wherefrom )

            #add cuts
            con = @build_constraint( model[:u][w] <= b.Phi + sum( b.rho[i] * ( model[:x][i] - x_val[i] ) for i ∈ 1:inst.nnodes ) )
            MOI.submit(model, MOI.UserCut(cb_data), con)
            res.nBENCUTS += 1
            res.nUCUTS   += 1
        end
    end # ∀ω ∈ Ω

    return
end


#****************************************************************
#**** IMP MODEL OUTER APPROXIMATION
#****************************************************************
# BUILD MODEL OA
function buildModel_OA( inst::instance, params::Dict, model::Model )
    print("build model...")

    #vars
    @variable( model, u[w=1:inst.nsc] >= 0 )                                            #scenario contribution
    @variable( model, x[i=1:inst.nnodes], Bin, base_name="x" )                          #seed sets
    @variable( model, xv[i=1:inst.nnodes, w=1:inst.nsc] >= 0, base_name="xv" )          #viewing variables
    @variable( model, 0 <= xf[i=1:inst.nnodes, w=1:inst.nsc] <= 1, base_name="xf" )     #forwarding variables

    #obj
    @objective( model, Max, (1/inst.nsc) * sum( u ) )

    #constraints
    @constraint( model, cardinality, sum( x[i] for i=1:inst.nnodes ) <= inst.kF )            # |F| <= k
    @constraint( model, forwarding[i=1:inst.nnodes, w=1:inst.nsc], sum( xf[j,w] for j in inst.lgi[i,w] ) >=  xv[i,w] )
    @constraint( model, viewing[i=1:inst.nnodes, w=1:inst.nsc], sum( x[j] for j in inst.Rin[i,w] ) >= xf[i,w] )
    @constraint( model, force_zero[i=1:inst.nnodes, w=1:inst.nsc], xv[i,w] <= length(inst.lgi[i,w]) * ( 1 - x[i] ) )

    # force variable values to zero for singletons
    for w in 1:inst.nsc
        for i in 1:inst.nnodes
            if inst.Singleton[i,w]          #if node i is a singleton in w
                @constraint( model, xv[i,w] == 0 )
                @constraint( model, xf[i,w] == 0 )
            end
        end
    end

    #initial cuts
    @constraint( model, init_cuts[w=1:inst.nsc], u[w] <= sum( (1 / inst.r[i] ) * xv[i,w] for i=1:inst.nnodes) )

    print("done! \n")
    return
end


# OUTER APPROXIMATION CALLBACK
function startOACallbackIMP( cb_data::CPLEX.CallbackContext, inst::instance, params::Dict, res::results, b::bfsdata, model::Model, wherefrom::String )
    #eps::Float64 = 0.01
    relcuttol::Float64 = params["relcuttol"]

    b.xv_val = annotateMatrix( callback_value.(Ref(cb_data), model[:xv]) )
    u_val = annotate( callback_value.(Ref(cb_data), model[:u]) )


    for w ∈ 1:inst.nsc
        clear!(b)
        for i in 1:inst.nnodes
            b.Phi += b.xv_val[i,w] / (b.xv_val[i,w] + inst.r[i])
        end


        if (u_val[w] / b.Phi) > relcuttol
            #add cuts
            con = @build_constraint( model[:u][w] <= b.Phi + sum( (inst.r[i] / (b.xv_val[i,w] + inst.r[i])^2) * (model[:xv][i,w] - b.xv_val[i,w] ) for i ∈ 1:inst.nnodes ) )
            if wherefrom == "LCB"
                MOI.submit(model, MOI.LazyConstraint(cb_data), con)
                res.nLCUTS   += 1
            end
            if wherefrom == "UCB"
                MOI.submit(model, MOI.UserCut(cb_data), con)
                res.nUCUTS   += 1
            end
            res.nBENCUTS += 1
        end
    end # ∀ω ∈ Ω

    return
end
