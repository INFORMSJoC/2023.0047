#****************************************************************
#**** MODELS CIMP GENERALIZED BENDERS DECOMPOSITION
#****************************************************************
function getInitialBendersCutsCIMP( inst::instance )
    ρ::Matrix{Float64}  = zeros( Float64, inst.nnodes, inst.nsc )                #cut coefficients
    α::Matrix{Float64}  = zeros( Float64, inst.nnodes, inst.nsc )                #the same alphas ∀ ω ∈ Ω

    #compute α[i,w] ∀ ω ∈ Ω
    for w=1:inst.nsc
        if inst.mtype == 5      #fractional benders
            for i=1:inst.nnodes
                α[i,w] = 1 / ( inst.rw[i,w] )
            end
        end
        if inst.mtype in [2,3]      #linear benders total/organic views
            for i=1:inst.nnodes
                α[i,w] = 1
            end
        end
    end

    #compute ρ[i,w]
    if inst.mtype in [2,3,5]
        for w=1:inst.nsc
            for i=1:inst.nnodes
                for k ∈ inst.Rout[i,w]
                    for j ∈ inst.lgo[k,w]
                        ρ[i,w] += α[j,w]
                    end
                end
                if i ∈ inst.SL
                    @assert( ρ[i,w] == 0, "ρ[$i,$w] != 0 in initial Benders cut")
                end
            end
        end
    end

    if inst.mtype == 1
        for w=1:inst.nsc
            for i=1:inst.nnodes
                for j ∈ inst.Rout[i,w]
                    ρ[i,w] += 1             #NOTE: β_j = 1
                end
                if i ∈ inst.SL
                    @assert( ρ[i,w] == 0, "ρ[$i,$w] != 0 in initial Benders cut")
                end
            end
        end
    end

    return ρ::Matrix{Float64}
end


# MODEL GENERALIZED BENDERS DECOMPOSITION
function buildModel_BENCIMP( inst::instance, params::Dict, model::Model )
    print("build model...")

    # get initial cut coefficients
    ρ::Matrix{Float64} = getInitialBendersCutsCIMP( inst )

    #vars
    @variable( model, u[w=1:inst.nsc] >= 0 )                           #scenario contribution
    @variable( model, y[i=1:inst.nnodes], Bin, base_name="y" )         #seed sets

    #obj
    @objective( model, Max, (1/inst.nsc) * sum( u ) )

    #constraints
    @constraint( model, sum( y[i] for i=1:inst.nnodes ) <= inst.kF )            # |F| <= k

    for w=1:inst.nsc
        @constraint( model, u[w] <= sum( ρ[i,w] * y[i] for i=1:inst.nnodes) )   # Benders Cut with x=0 ∀i ∈ V
    end

    for i ∈ inst.SL
        @constraint( model, y[i] == 0 )
    end

    print("done! \n")
    return
end


#****************************************************************
#**** CIMP CALLBACKS
#****************************************************************
#**** LAZY 

function computeSubgradientCIMP!( inst::instance, res::results, b::bfsdata, S::Vector{Tnodes}, w::Int, wherefrom::String )
    #compute α_i / φ_i
    if inst.mtype == 5                                  #GBD
        for i ∈ 1:inst.nnodes
            if !inst.Singleton[i,w]                     #dual vars do not exist for singletons
                if ( b.yv[i] >= inst.indeg[i,w] )
                    b.phi[i] = (inst.rw[i,w]) / (inst.indeg[i,w] + inst.rw[i,w])^2
                else
                    b.alpha[i] = (inst.rw[i,w]) / (b.yv[i] + inst.rw[i,w])^2
                end
            end
        end
        for i ∈ S
            b.phi[i] = 1 / (inst.rw[i,w])
            b.alpha[i] = 0.0
        end
        for i ∈ inst.SL
            b.alpha[i] = 0.0
            b.phi[i] = 0.0
        end
    end

    if inst.mtype in [2,3]                              #total views / organic reach
        for i ∈ 1:inst.nnodes
            if !inst.Singleton[i,w]                     #dual vars do not exist for singletons
                if ( b.yv[i] >= inst.indeg[i,w] )
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
        for i ∈ inst.SL
            b.alpha[i] = 0.0
            b.phi[i] = 0.0
        end
    end

    #compute ρ_i
    if inst.mtype in [2,3,5]
        for i ∈ 1:inst.nnodes
            if !inst.Singleton[i,w]                     #dual vars do not exist for singletons
                for k ∈ inst.Rout[i,w]
                    if !b.yfb[k]                              #β_k = 0 ∀ forwarders
                        for j ∈ inst.lgo[k,w]
                            b.rho[i] += b.alpha[j]
                        end
                    end
                end
                b.rho[i] -= inst.indeg[i,w] * b.phi[i]
            end
        end
    
        for i in inst.SL
            b.rho[i] = 0.0
        end
    end

    if inst.mtype == 1          #forwarders
        if wherefrom == "LCB"
            for i ∈ 1:inst.nnodes
                if ( b.yfb[i] == 0 )
                    b.alpha[i] = 1
                end
            end
        end
        if wherefrom == "UCB"
            for i ∈ 1:inst.nnodes
                if ( b.yf[i] < 1 )
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


function computeLazyObjectiveCIMP!( inst::instance, res::results, b::bfsdata, S::Vector{Tnodes}, w::Int )
    #compute yf, yv
    for i ∈ S
        for k ∈ inst.Rout[i,w]
            if !b.yfb[k]                     #forwarding not marked yet
                b.yfb[k] = true
                for j ∈ inst.lgo[k,w]
                    b.yv[j] += 1.0
                end
            end
        end
    end
    for i ∈ S
        b.yv[i] = 0.0
    end
    for i in inst.SL
        b.yv[i] = 0.0
    end

    #adjust viewings for organic case
    if inst.mtype == 3
        for i in 1:inst.nnodes
            if b.yv[i] > 0
                b.yv[i] = 1.0
            end
        end
    end

    #compute solution of Benders subproblem
    if inst.mtype == 5      #fractional benders
        for i ∈ 1:inst.nnodes
            b.Phi  += b.yv[i] / (b.yv[i] + inst.rw[i,w])
        end
    end
    if inst.mtype in [2,3]      #linear benders total/organic views
        for i ∈ 1:inst.nnodes
            b.Phi  += b.yv[i]
        end
    end
    if inst.mtype == 1      #Kempe
        for i ∈ 1:inst.nnodes
            b.Phi  += b.yfb[i]
        end
    end

    return
end


#LAZY CUTS GBD
function startLazyBendersCallbackCIMP( cb_data::CPLEX.CallbackContext, inst::instance, params::Dict, res::results, b::bfsdata, model::Model )
    relcuttol::Float64 = params["relcuttol"]
    wherefrom::String = "LCB"

    #retrieve vals
    y_val = annotate( callback_value.(Ref(cb_data), model[:y]) )
    u_val = annotate( callback_value.(Ref(cb_data), model[:u]) )

    # set current seed set
    S = Tnodes[]
    sizehint!(S, inst.kF)
    for i ∈ eachindex( y_val )
        if y_val[i] > 0.5
            push!(S, i)
        end
    end

    for w ∈ 1:inst.nsc
        clear!(b)
        computeLazyObjectiveCIMP!( inst, res, b, S, w )

        if (u_val[w] / b.Phi) > relcuttol
            computeSubgradientCIMP!( inst, res, b, S, w, wherefrom )

            #add cuts
            con = @build_constraint( model[:u][w] <= b.Phi + sum( b.rho[i] * model[:y][i] for i ∈ 1:inst.nnodes ) - sum( b.rho[i] for i ∈ S ) )
            MOI.submit(model, MOI.LazyConstraint(cb_data), con)
            res.nBENCUTS += 1
            res.nLCUTS   += 1
        end
    end # ∀ω ∈ Ω

    return
end
#***END LAZY 

#***USER 
function computeUserObjectiveCIMP!( inst::instance, res::results, b::bfsdata, y_val::Vector{Float64}, w::Int )
    tempsum::Float64 = 0.0

    #compute yf
    for i ∈ 1:inst.nnodes
        if !inst.Singleton[i,w]                         #all values are zero for singletons
            tempsum = 0.0
            for k ∈ inst.Rin[i,w]
                tempsum += y_val[k]
                tempsum >= 1.0 ? break : nothing
            end
            b.yf[i] = min( 1.0, tempsum )
            if b.yf[i] == 1.0                 
                b.yfb[i] = true
            end
        end
    end

    for i in inst.SL
        #b.xv[i]  = 0
        b.yf[i]  = 0
        b.yfb[i] = false
    end

    #compute yv, π, Φ(y)
    if inst.mtype == 5      #fractional benders
        for i ∈ 1:inst.nnodes
            if !inst.Singleton[i,w]                         #all values are zero for singletons
                tempsum = 0.0
                for j ∈ inst.lgi[i,w]
                    tempsum += b.yf[j]
                end
                b.yv[i] = min( inst.indeg[i,w] * (1-y_val[i]), tempsum )
                if i in inst.SL
                    b.yv[i] = 0
                end
                b.Phi  += b.yv[i] / ( b.yv[i] + inst.rw[i,w])
            end
        end
    end

    if inst.mtype == 2      # total views
        for i ∈ 1:inst.nnodes
            if !inst.Singleton[i,w]                         #all values are zero for singletons
                tempsum = 0.0
                for j ∈ inst.lgi[i,w]
                    tempsum += b.yf[j]
                end
                b.yv[i] = min( inst.indeg[i,w] * (1-y_val[i]), tempsum )
                if i in inst.SL
                    b.yv[i] = 0
                end
                b.Phi  += b.yv[i]
            end
        end
    end

    if inst.mtype == 3      # organic reach
        for i ∈ 1:inst.nnodes
            if !inst.Singleton[i,w]                         #all values are zero for singletons
                tempsum = 0.0
                for j ∈ inst.lgi[i,w]
                    tempsum += b.yf[j]
                end
                b.yv[i] = min( inst.indeg[i,w] * (1-y_val[i]), tempsum )
                if i in inst.SL
                    b.yv[i] = 0
                end
                if b.yv[i] > 1
                    b.yv[i] = 1.0
                end
                b.Phi  += b.yv[i]
            end
        end
    end

    if inst.mtype == 1                  # forwards
        for i ∈ 1:inst.nnodes
            b.Phi  += b.yf[i]
        end
    end

    return
end


#USER BENDERS CALLBACK CIMP
function startUserBendersCallbackCIMP( cb_data::CPLEX.CallbackContext, inst::instance, params::Dict, res::results, b::bfsdata, model::Model )
    eps::Float64 = 0.001
    relcuttol::Float64 = deepcopy( params["relcuttolU"] )
    wherefrom::String = "UCB"

    # get current values
    y_val = annotate( callback_value.(Ref(cb_data), model[:y]) )
    u_val = annotate( callback_value.(Ref(cb_data), model[:u]) )


    # set current seed set
    S::Vector{Tnodes} = Tnodes[]
    sizehint!(S, inst.kF)
    for i ∈ eachindex( y_val )
        if y_val[i] > 1-eps
            y_val[i] = 1.0
            push!(S, i)
        end
        if y_val[i] < eps
            y_val[i] = 0.0
        end
    end

    #compute and add cuts
    for w ∈ 1:inst.nsc
        clear!(b)
        computeUserObjectiveCIMP!( inst, res, b, y_val, w )

        if (u_val[w] / b.Phi) > relcuttol
            computeSubgradientCIMP!( inst, res, b, S, w, wherefrom )

            #add cuts
            con = @build_constraint( model[:u][w] <= b.Phi + sum( b.rho[i] * ( model[:y][i] - y_val[i] ) for i ∈ 1:inst.nnodes ) )
            MOI.submit(model, MOI.UserCut(cb_data), con)
            res.nBENCUTS += 1
            res.nUCUTS   += 1
        end
    end # ∀ω ∈ Ω

    return
end
#***END USER VERSION



#****************************************************************
#**** CIMP MODEL OUTER APPROXIMATION
#****************************************************************
# BUILD MODEL
function buildModel_OACIMP( inst::instance, params::Dict, model::Model )
    print("build model...")

    #vars
    @variable( model, u[w=1:inst.nsc] >= 0 )                                        #scenario contribution
    @variable( model, y[i=1:inst.nnodes], Bin, base_name="y" )                      #seed sets
    @variable( model, yv[i=1:inst.nnodes, w=1:inst.nsc] >= 0, base_name="yv" )      #viewing variables
    @variable( model, 0 <= yf[i=1:inst.nnodes, w=1:inst.nsc] <= 1, base_name="yf" ) #forwarding variables      

    #obj
    @objective( model, Max, (1/inst.nsc) * sum( u ) )

    #constraints
    @constraint( model, cardinality, sum( y[i] for i=1:inst.nnodes ) == inst.kF )            # |F| <= k
    @constraint( model, forwarding[i=1:inst.nnodes, w=1:inst.nsc], sum( yf[j,w] for j in inst.lgi[i,w] ) >=  yv[i,w] )
    @constraint( model, viewing[i=1:inst.nnodes, w=1:inst.nsc], sum( y[j] for j in inst.Rin[i,w] ) >= yf[i,w] )
    @constraint( model, force_zero[i=1:inst.nnodes, w=1:inst.nsc], yv[i,w] <= length(inst.lgi[i,w]) * ( 1 - y[i] ) )

    # set variables to zero for singeltons and leader seed nodes
    for w in 1:inst.nsc
        for i in 1:inst.nnodes
            if inst.Singleton[i,w]          #if node i is a singleton in w
                @constraint( model, yv[i,w] == 0 )
                @constraint( model, yf[i,w] == 0 )
            end
        end
    end

    for w in 1:inst.nsc
        for i in inst.SL
            @constraint( model, yv[i,w] == 0 )
            @constraint( model, yf[i,w] == 0 )
        end
    end

    #initial cuts
    @constraint( model, init_cuts[w=1:inst.nsc], u[w] <= sum( (1 / inst.rw[i,w] ) * yv[i,w] for i=1:inst.nnodes) ) 

    print("done! \n")
    return
end


# OUTER APPROXIMATION CALLBACK
function startOACallbackCIMP_forallRHS( cb_data::CPLEX.CallbackContext, inst::instance, params::Dict, res::results, b::bfsdata, model::Model, wherefrom::String )
    relcuttol::Float64 = params["relcuttol"]

    b.yv_val = annotateMatrix( callback_value.(Ref(cb_data), model[:yv]) )
    u_val = annotate( callback_value.(Ref(cb_data), model[:u]) )

    for w ∈ 1:inst.nsc
        clear!(b)
        for i in 1:inst.nnodes
            b.Phi += b.yv_val[i,w] / (b.yv_val[i,w] + inst.rw[i,w])
        end

        #add violated cuts
        if (u_val[w] / b.Phi) > relcuttol
            con = @build_constraint( model[:u][w] <= b.Phi + sum( (inst.rw[i,w] / (b.yv_val[i,w] + inst.rw[i,w])^2) * (model[:yv][i,w] - b.yv_val[i,w] ) for i ∈ 1:inst.nnodes ) )
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

