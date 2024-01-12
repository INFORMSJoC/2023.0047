#****************************************************************
#**** FUNCTIONS RELATED TO GRAPHS (AND MORE)
#****************************************************************

# CREATE LIVE-GRAPHS (OR SCENARIO GRAPHS)
function createLiveGraphs!( inst::instance )
    #create empty array of sets
    inst.lgo   = [Tnodes[] for i=1:inst.nnodes, j=1:inst.nsc]
    inst.lgi   = [Tnodes[] for i=1:inst.nnodes, j=1:inst.nsc]
    inst.lgfo  = [Tnodes[] for i=1:inst.nnodes, j=1:inst.nsc]
    inst.lgfi  = [Tnodes[] for i=1:inst.nnodes, j=1:inst.nsc]
    inst.indeg = zeros( Tnodes, inst.nnodes, inst.nsc )
    inst.nlivearcs = zeros(Int, inst.nsc)
    i::Tnodes = 0
    j::Tnodes = 0
    pf::Tprob = 0.0
    pv::Tprob = 0.0
    r::Tprob  = 0.0

    aid::Int64 = 0              # dummy arc id
    for w in 1:inst.nsc
        aid = 0
        for a ∈ inst.edge
            aid +=1
            i, j, pf, pv = a
            if( pf == pv )                              #empirical pf > params["probV"]
                if( rand( Bernoulli( pf ) ) == 1 )
                    inst.indeg[j, w] += 1
                    push!( inst.lgo[i, w], j )
                    push!( inst.lgi[j, w], i )
                    push!( inst.lgfo[i, w], j )
                    push!( inst.lgfi[j, w], i )
                    inst.nlivearcs[w] += 1
                end
            else
                r = rand(Tprob)
                if(r <= pf )
                    inst.indeg[j, w] += 1
                    push!( inst.lgo[i, w], j )
                    push!( inst.lgi[j, w], i )
                    push!( inst.lgfo[i, w], j )
                    push!( inst.lgfi[j, w], i )
                    inst.nlivearcs[w] += 1
                elseif( pf < r <= pv )
                    inst.indeg[j, w] += 1
                    push!( inst.lgo[i, w], j )
                    push!( inst.lgi[j, w], i )
                    inst.nlivearcs[w] += 1
                else
                    # do nothing
                end
            end
        end
        @assert( aid == inst.narcs, "aid != inst.narcs" )
        mb = get_mem_use()
        inst.memok = memOK(mb, inst.memlimit)
        if !inst.memok
            return
        end
    end # ∀ω ∈ Ω

    if inst.mtype in [1,3] # forwarding / organic views
        inst.indeg = ones( Tnodes, inst.nnodes, inst.nsc )
    end

    return nothing
end # createLiveGraphs

# return number elements in live graph (note: not equal to number of live-arcs)
function getNElementsInLiveGraph( inst::instance )
    nelements::Int64 = 0
    for w=1:inst.nsc
        for i=1:inst.nnodes
            nelements += length( inst.lgo[i,w] )
        end
    end
    return nelements::Int64
end


# COMPUTE RESISTANCE VALUES  r_i, ∀ i ∈ V
function getResistanceValues!( inst::instance, params::Dict )
    inst.r              = zeros( Float64, inst.nnodes )
    multiplier::Float64 = (1- params["resistanceHurdle"]) / params["resistanceHurdle"]

    #determine which nodes are resistant
    is_resistant      = falses(inst.nnodes)          #indicates if a node is resistant or not
    resistant_ranks   = Tuple{Int,Float64}[]         #criteria which nodes are resistant and which not
    n_resistant_nodes = round(Int, inst.nnodes * params["fractionOfResistantNodes"])
    if params["fractionOfResistantNodes"] > 0
        for i in 1:inst.nnodes
            push!(resistant_ranks, (i, rand()))
        end
        sort!(resistant_ranks, by = x->(x[2]) )

        for i in 1:n_resistant_nodes
            nodeID, val = resistant_ranks[i]
            is_resistant[nodeID] = true
        end
    end

    # compute resistance values
    for i=1:inst.nnodes
        if indegree(inst.g, i) != 0
            if is_resistant[i]              #node is very resistant
                inst.r[i] = round(Int, multiplier * indegree(inst.g, i) )
            else
                inst.r[i] = rand( 1 : 1 : (multiplier * indegree(inst.g, i))/2 )
            end
        else
            inst.r[i] = 1
        end
    end

    return nothing
end

# COMPUTE SCENARIO DEPENDENT RESISTANCE VALUES
function compute_rw!( inst::instance, b::bfsdata )
    inst.rw = zeros( Float64, inst.nnodes, inst.nsc )

    for w in 1:inst.nsc
        clear!(b)

        for i in 1:inst.nnodes
            inst.rw[i,w] = inst.r[i]
        end
        # set seed sets to current level
        for i ∈ inst.SL
            b.bvisited[i] = true
            push!(b.cur_level_L, i)
        end

        #propagation
        while !( isempty(b.cur_level_L) )
            #propagate from L
            @inbounds for v in b.cur_level_L
                #increment outneigbors
                @inbounds for i in inst.lgo[v,w] #arcs (v,i), i.e., outneighbors
                    inst.rw[i,w] += 1 * inst.leaderUtility
                end
                #propagate along forwarding arcs
                @inbounds for i in inst.lgfo[v,w] #arcs (v,i), i.e., outneighbors
                        if( b.bvisited[i] == false )
                            push!(b.next_level_L, i)
                            b.bvisited[i] = true
                        end
                end
            end
            empty!(b.cur_level_L)
            b.cur_level_L, b.next_level_L = b.next_level_L, b.cur_level_L
        end

        for i in inst.SL
            inst.rw[i,w] = inst.r[i]
        end
    end

    return nothing
end


#****************************************************************
#**** (REVERSE) ACTIVATION SETS
#****************************************************************
# COMPUTE REVERSE ACTIVATION SET FOR A SCENARIO AND NODE (IMP)
function propagateValidActivatorsIMP!( inst::instance, b::bfsdata, w::Int, j::Int )
    #clear bfs data
    clear!(b)

    # set seed sets to current level
    b.bvisited[j] = true
    push!(b.cur_level_F, j)
    push!( inst.Rin[j,w], j )

    #propagation
    while ( !isempty(b.cur_level_F) )
        @inbounds for v in b.cur_level_F
            @inbounds for i in inst.lgfi[v,w] #arcs (i,v), i.e., inneighbors
                if( b.bvisited[i] == false )
                    @views if( i < j ) # add nodes from R_i to R_j (since R_i was already comuted)
                                for k ∈ inst.Rin[i,w]
                                    if( !b.bvisited[k] )
                                        b.bvisited[k] = true
                                        push!( inst.Rin[j,w], k )
                                    end
                                end
                            else        #not computed R^w_i; so just add node i
                                push!(b.next_level_F, i)
                                b.bvisited[i] = true
                                push!( inst.Rin[j,w], i )
                            end #endif @views
                end #endif b.visited == false
            end # for i
        end # for v
        empty!(b.cur_level_F)
        b.cur_level_F, b.next_level_F = b.next_level_F, b.cur_level_F
    end #while


    #set nodes that are reach node j
    for k ∈ inst.Rin[j,w]
        push!( inst.Rout[k,w], j)
    end
    return nothing
end


# COMPUTE VALID ACTIVATORS IMP (i.e., reverse activation set  A^-_ω(i), ∀ ω ∈ Ω)
function getValidActivatorsIMP!( inst::instance, b::bfsdata )
    inst.Rin  = [Tnodes[] for i=1:inst.nnodes, j=1:inst.nsc]
    inst.Rout = [Tnodes[] for i=1:inst.nnodes, j=1:inst.nsc]
    inst.Singleton   = falses( inst.nnodes, inst.nsc )
    inst.isSingleton = zeros( Bool, inst.nnodes )

    #singleton check
    has_outneighbor_f::Bool = false
    has_inneighbor_f::Bool  = false
    has_outneighbor::Bool   = false
    has_inneighbor::Bool    = false

    for w=1:inst.nsc
        for j=1:inst.nnodes
            isempty(inst.lgo[j,w]) ? nothing : has_outneighbor = true
            isempty(inst.lgi[j,w]) ? nothing : has_inneighbor  = true
            isempty(inst.lgfo[j,w]) ? nothing : has_outneighbor_f = true
            isempty(inst.lgfi[j,w]) ? nothing : has_inneighbor_f  = true
            if(!has_inneighbor && !has_outneighbor)
                inst.Singleton[j,w] = 1
                inst.isSingleton[j] = 1
            end
            if( has_inneighbor_f )
                propagateValidActivatorsIMP!( inst, b, w, j )
            else
                push!( inst.Rin[j,w], j )
                push!( inst.Rout[j,w], j )
            end
            has_outneighbor_f = false
            has_inneighbor_f  = false
            has_outneighbor   = false
            has_inneighbor    = false
        end
        #check memory
        mb = get_mem_use()
        inst.memok = memOK(mb, inst.memlimit)
        if !inst.memok
            return
        end
    end

    return nothing
end


# PROPAGATE A^-_ω(j) CIMP
function propagateValidActivatorsCIMP!( inst::instance, b::bfsdata, w::Int, j::Int )
    #clear bfs data
    clear!(b)

    # set seed sets to current level
    b.bvisited[j] = true
    push!(b.cur_level_F, j)
    push!( inst.Rin[j,w], j )

    #propagation
    while ( !isempty(b.cur_level_F) )
        @inbounds for v in b.cur_level_F
            @inbounds for i in inst.lgfi[v,w] #arcs (i,v), i.e., inneighbors
                if( b.bvisited[i] == false )
                    @views if( i < j ) # add nodes from A_i to A_j (since A_i was already comuted)
                                for k ∈ inst.Rin[i,w]
                                    if( !b.bvisited[k] )
                                        b.bvisited[k] = true
                                        push!( inst.Rin[j,w], k )
                                    end
                                end
                            else        #not computed A^w_i; so just add node i
                                b.bvisited[i] = true
                                if !inst.inSL[i]
                                    push!(b.next_level_F, i)
                                    push!( inst.Rin[j,w], i )
                                end
                            end #endif @views
                end #endif b.visited == false
            end # for i
        end # for v
        empty!(b.cur_level_F)
        b.cur_level_F, b.next_level_F = b.next_level_F, b.cur_level_F
    end #while


    #set nodes that are reach node j
    for k ∈ inst.Rin[j,w]
        push!( inst.Rout[k,w], j)
    end
    return nothing
end


# CIMP: COMPUTE VALID ACTIVATORS A^-_ω(j) ∀ ω ∈ Ω
function getValidActivatorsCIMP!( inst::instance, b::bfsdata )
    inst.Rin  = [Tnodes[] for i=1:inst.nnodes, j=1:inst.nsc]
    inst.Rout = [Tnodes[] for i=1:inst.nnodes, j=1:inst.nsc]
    inst.Singleton   = falses( inst.nnodes, inst.nsc )
    inst.isSingleton = zeros( Bool, inst.nnodes )

    #singleton check
    has_outneighbor_f::Bool = false
    has_inneighbor_f::Bool  = false
    has_outneighbor::Bool   = false
    has_inneighbor::Bool    = false

    for w=1:inst.nsc
        for j=1:inst.nnodes
            if !inst.inSL[j]
                isempty(inst.lgo[j,w]) ? nothing : has_outneighbor = true
                isempty(inst.lgi[j,w]) ? nothing : has_inneighbor  = true
                isempty(inst.lgfo[j,w]) ? nothing : has_outneighbor_f = true
                isempty(inst.lgfi[j,w]) ? nothing : has_inneighbor_f  = true
                if(!has_inneighbor && !has_outneighbor)
                    inst.Singleton[j,w] = 1
                    inst.isSingleton[j] = 1
                end
                if( has_inneighbor_f )
                    propagateValidActivatorsCIMP!( inst, b, w, j )
                else
                    push!( inst.Rin[j,w], j )
                    push!( inst.Rout[j,w], j )
                end
                has_outneighbor_f = false
                has_inneighbor_f  = false
                has_outneighbor   = false
                has_inneighbor    = false
            end
        end
        #check memory
        mb = get_mem_use()
        inst.memok = memOK(mb, inst.memlimit)
        if !inst.memok
            return
        end
    end

    return nothing
end


function getNElementsInReachabilityset( inst::instance )
    nelements::Int64 = 0
    for w=1:inst.nsc
        for i=1:inst.nnodes
            nelements += length( inst.Rin[i,w] )
        end
    end
    return nelements::Int64
end


#****************************************************************
#**** COMPUTE OBJECTIVE FUNCTIONS FROM SET (small set of scenarios Ω' ⊂ Ω)
#****************************************************************
# IMP
function computeObjectiveValuesFromSetIMP( inst::instance, b::bfsdata, S::Array{Int,1} )
    obj_val::Float64 = 0.0
    u = zeros(Float64, inst.nsc)    #contribution of each scenario (μ^ω)

    for w = 1:inst.nsc
        clear!(b)
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

        #adjust viewings for organic case
        if inst.mtype == 3
            for i in 1:inst.nnodes
                if b.xv[i] > 1.0
                    b.xv[i] = 1.0
                end
            end
        end

        #compute solution of Benders subproblem / alphas
        if inst.mtype in [4,5]
            for i ∈ 1:inst.nnodes
                u[w]  += b.xv[i] / (b.xv[i] + inst.r[i])
            end
        end
        if inst.mtype in [2,3]
            for i ∈ 1:inst.nnodes
                u[w]  += b.xv[i]
            end
        end
        if inst.mtype == 1
            for i ∈ 1:inst.nnodes
                u[w]  += b.xfb[i]
            end
        end

    end
    obj_val = sum( u ) * (1/inst.nsc)

    return obj_val::Float64
end

# CIMP
function computeObjectiveValuesFromSetCIMP( inst::instance, b::bfsdata, S::Array{Int,1} )
    obj_val::Float64 = 0.0
    u = zeros(Float64, inst.nsc)    #contribution of each scenario (μ^ω)

    for w = 1:inst.nsc
        clear!(b)
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

        #adjust viewings for organic views
        if inst.mtype == 3
            for i in 1:inst.nnodes
                if b.yv[i] > 1
                    b.yv[i] = 1.0
                end
            end
        end

        #compute solution of Benders subproblem / alphas
        if inst.mtype == 5      #fractional benders
            for i ∈ 1:inst.nnodes
                u[w]  += b.yv[i] / (b.yv[i] + inst.rw[i,w])
            end
        end
        if inst.mtype in [2,3]      #linear benders
            for i ∈ 1:inst.nnodes
                u[w]  += b.yv[i]
            end
        end
        if inst.mtype == 1          #Kempe
            for i ∈ 1:inst.nnodes
                u[w]  += b.yfb[i]
            end
        end
    end
    obj_val = sum( u ) * (1/inst.nsc)

    return obj_val::Float64
end


#****************************************************************
#**** EVALUATE SEED SETS (on large set of scenarios Ω'' ⊂ \Omega)
#****************************************************************
# IMP
function computeViewingValuesEvalIMP!( inst::instance, b::bfsdata, SL::Array{Int,1}, SF::Array{Int,1} )
    # set seed sets to current level
    for i ∈ SF
        b.bvisited[i] = true
        push!(b.cur_level_F, i)
    end

    #create live graphs
    while !( isempty(b.cur_level_F) )
        #propagate from F
        @inbounds for i in b.cur_level_F
                        for t in inst.gout[i]
                            j, pf, pv = t[1], t[2], t[3]
                            if( pf == pv )                              #empirical pf > params["probV"]
                                if( rand( Bernoulli( pf ) ) == 1 )
                                    push!( inst.lgo_tmp[i], j )
                                    push!( inst.lgfo_tmp[i], j )
                                    if !b.bvisited[j]
                                        b.bvisited[j] = true
                                        push!(b.next_level_F, j)
                                    end
                                end
                            else
                                r = rand(Tprob)
                                if(r <= pf )
                                    push!( inst.lgo_tmp[i], j )
                                    push!( inst.lgfo_tmp[i], j )
                                    if !b.bvisited[j]
                                        b.bvisited[j] = true
                                        push!(b.next_level_F, j)
                                    end
                                elseif( pf < r <= pv )
                                    push!( inst.lgo_tmp[i], j )
                                else
                                    # do nothing
                                end
                            end
                        end
                    end
                empty!(b.cur_level_F)
                b.cur_level_F, b.next_level_F = b.next_level_F, b.cur_level_F
    end

    #propagate
    fill!(b.bvisited, zero(Bool))
    for i ∈ SF
        b.bvisited[i] = true
        push!(b.cur_level_F, i)
    end

    while !( isempty(b.cur_level_F) )
        #propagate from F
        @inbounds for v in b.cur_level_F
                      b.yf[v] = 1.0
            @inbounds for k in inst.lgfo_tmp[v]
                        if !b.bvisited[k]
                            b.bvisited[k] = true
                            push!(b.next_level_F, k)
                        end
                      end
                      for k in inst.lgo_tmp[v]
                          b.yv[k] += 1.0
                      end
                   end
        empty!(b.cur_level_F)
        b.cur_level_F, b.next_level_F = b.next_level_F, b.cur_level_F
    end

    for i ∈ SF
        b.yv[i] = 0
    end

    for i=1:inst.nnodes
        empty!(inst.lgo_tmp[i])
        empty!(inst.lgfo_tmp[i])
    end

    return
end


# CIMP
function computeViewingValuesEvalCIMP!( inst::instance, b::bfsdata, SL::Array{Int,1}, SF::Array{Int,1} )
    # set seed sets to current level
    for i ∈ SF
        b.bvisited[i] = true
        push!(b.cur_level_F, i)
    end

    for i ∈ SL
        b.bvisited[i] = true
        push!(b.cur_level_F, i)             #NOTE: for live graph creation both leader and follower are pushed in one BFS level
    end


    #create live graphs
    while !( isempty(b.cur_level_F) )
        #propagate from F
        @inbounds for i in b.cur_level_F
                        for t in inst.gout[i]
                            j, pf, pv = t[1], t[2], t[3]
                            if( pf == pv )                              #empirical pf > params["probV"]
                                if( rand( Bernoulli( pf ) ) == 1 )
                                    push!( inst.lgo_tmp[i], j )
                                    push!( inst.lgfo_tmp[i], j )
                                    if !b.bvisited[j]
                                        b.bvisited[j] = true
                                        push!(b.next_level_F, j)
                                    end
                                end
                            else
                                r = rand(Tprob)
                                if(r <= pf )
                                    push!( inst.lgo_tmp[i], j )
                                    push!( inst.lgfo_tmp[i], j )
                                    if !b.bvisited[j]
                                        b.bvisited[j] = true
                                        push!(b.next_level_F, j)
                                    end
                                elseif( pf < r <= pv )
                                    push!( inst.lgo_tmp[i], j )
                                else
                                    # do nothing
                                end
                            end
                        end
                    end
                empty!(b.cur_level_F)
                b.cur_level_F, b.next_level_F = b.next_level_F, b.cur_level_F
    end

    #propagate from L
    fill!(b.bvisited, zero(Bool))
    for i ∈ SL
        b.bvisited[i] = true
        push!(b.cur_level_L, i)
    end

    while !( isempty(b.cur_level_L) )
        #propagate from L
        @inbounds for v in b.cur_level_L
                      b.xf[v] = 1.0
            @inbounds for k in inst.lgfo_tmp[v]
                        if !b.bvisited[k]
                            b.bvisited[k] = true
                            push!(b.next_level_L, k)
                        end
                      end
                      for k in inst.lgo_tmp[v]
                          b.xv[k] += 1.0
                      end
                   end
        empty!(b.cur_level_L)
        b.cur_level_L, b.next_level_L = b.next_level_L, b.cur_level_L
    end

    #propagate from F
    fill!(b.bvisited, zero(Bool))
    for i ∈ SF
        b.bvisited[i] = true
        push!(b.cur_level_F, i)
    end

    while !( isempty(b.cur_level_F) )
        @inbounds for v in b.cur_level_F
                      b.yf[v] = 1.0
            @inbounds for k in inst.lgfo_tmp[v]
                        if !b.bvisited[k]
                            b.bvisited[k] = true
                            if !inst.inSL[k]
                                push!(b.next_level_F, k)
                            end
                        end
                      end
                      for k in inst.lgo_tmp[v]
                          b.yv[k] += 1.0
                      end
                   end
        empty!(b.cur_level_F)
        b.cur_level_F, b.next_level_F = b.next_level_F, b.cur_level_F
    end

    for i ∈ SF
        b.xv[i] = 0
        b.yv[i] = 0
    end

    for i ∈ SL
        b.xv[i] = 0
        b.yv[i] = 0
    end

    for i=1:inst.nnodes
        empty!(inst.lgo_tmp[i])
        empty!(inst.lgfo_tmp[i])
    end

    return
end


# EVALUATE SEED SETS (on large set of scenarios Ω')
# NOTE: random seed is set before calling this function
function evaluateObjectiveFunction( inst::instance, b::bfsdata, SL::Array{Int,1}, SF::Array{Int,1}, wherefrom::String="OTHER_METHOD" )
    obj_val_F::Float64 = 0.0
    obj_val_L::Float64 = 0.0
    u::Float64 = 0.0
    uL::Float64 = 0.0
    u_vals = zeros(Float64, inst.Nsc)
    inst.lgo_tmp  = [Tnodes[] for i=1:inst.nnodes]
    inst.lgfo_tmp = [Tnodes[] for i=1:inst.nnodes]
    if wherefrom == "BENDERS"
        inst.nTotalViewsF       = zeros(Float64, inst.nnodes)
        inst.nOrganicViewsF     = zeros(Float64, inst.nnodes)
        inst.nForwardsF         = zeros(Float64, inst.nnodes)
        inst.nTotalViewsL       = zeros(Float64, inst.nnodes)
        inst.nOrganicViewsL     = zeros(Float64, inst.nnodes)
        inst.nForwardsL         = zeros(Float64, inst.nnodes)
    end

    for w = 1:inst.Nsc
        clear!(b)
        if inst.kL == 0
            computeViewingValuesEvalIMP!( inst, b, SL, SF )
        else
            computeViewingValuesEvalCIMP!( inst, b, SL, SF )
        end

        #adjust viewings for organic case
        if inst.mtype == 3
            for i in 1:inst.nnodes
                b.xv[i] > 1.0 ? b.xv[i] = 1.0 : nothing
                b.yv[i] > 1.0 ? b.yv[i] = 1.0 : nothing
            end
        end

        if wherefrom == "OTHER_METHOD"
            if inst.mtype in [4,5]
                for i in 1:inst.nnodes
                    u  += b.yv[i] / (inst.leaderUtility * b.xv[i] + b.yv[i] + inst.r[i])          #CAUTION: yv is also used for IMP here, thus xv=0
                    uL += (inst.leaderUtility * b.xv[i]) / (inst.leaderUtility * b.xv[i] + b.yv[i] + inst.r[i])          #CAUTION: yv is also used for IMP here, thus this is zero in the IMP
                    u_vals[w] += b.yv[i] / (b.xv[i] + b.yv[i] + inst.r[i])   #CAUTION: yv is also used for IMP here
                end
            end
            if inst.mtype in [2,3]
                for i in 1:inst.nnodes
                    u  += b.yv[i]               #CAUTION: yv is also used for IMP here, thus xv=0
                    uL += b.xv[i]               #CAUTION: yv is also used for IMP here, thus this is zero in the IMP
                    u_vals[w] += b.yv[i]        #CAUTION: yv is also used for IMP here
                end
            end
            if inst.mtype == 1
                for i in 1:inst.nnodes
                    u  += b.yf[i]               #CAUTION: yv is also used for IMP here, thus xv=0
                    uL += b.xf[i]               #CAUTION: yv is also used for IMP here, thus this is zero in the IMP
                    u_vals[w] += b.yf[i]        #CAUTION: yv is also used for IMP here
                end
            end
        end
        if wherefrom == "BENDERS"
            if inst.mtype in [4,5]
                for i in 1:inst.nnodes
                    u  += b.yv[i] / (inst.leaderUtility * b.xv[i] + b.yv[i] + inst.r[i])          #CAUTION: yv is also used for IMP here, thus xv=0
                    uL += (inst.leaderUtility * b.xv[i]) / (inst.leaderUtility * b.xv[i] + b.yv[i] + inst.r[i])          #CAUTION: yv is also used for IMP here, thus this is zero in the IMP
                    u_vals[w] += b.yv[i] / (b.xv[i] + b.yv[i] + inst.r[i])   #CAUTION: yv is also used for IMP here
                    inst.nTotalViewsF[i]       += b.yv[i]
                    inst.nOrganicViewsF[i]     += min(b.yv[i], 1)
                    inst.nForwardsF[i]         += b.yf[i]
                    inst.nTotalViewsL[i]       += b.xv[i]
                    inst.nOrganicViewsL[i]     += min(b.xv[i], 1)
                    inst.nForwardsL[i]         += b.xf[i]
                end
            end
            if inst.mtype in [2,3]
                for i in 1:inst.nnodes
                    u  += b.yv[i]               #CAUTION: yv is also used for IMP here, thus xv=0
                    uL += b.xv[i]               #CAUTION: yv is also used for IMP here, thus this is zero in the IMP
                    u_vals[w] += b.yv[i]        #CAUTION: yv is also used for IMP here
                    inst.nTotalViewsF[i]       += b.yv[i]
                    inst.nOrganicViewsF[i]     += min(b.yv[i], 1)
                    inst.nForwardsF[i]         += b.yf[i]
                    inst.nTotalViewsL[i]       += b.xv[i]
                    inst.nOrganicViewsL[i]     += min(b.xv[i], 1)
                    inst.nForwardsL[i]         += b.xf[i]
                end
            end
        end
    end

    obj_val_F = u * (1/inst.Nsc)
    obj_val_L = uL * (1/inst.Nsc)

    if wherefrom == "OTHER_METHOD"
        return obj_val_F::Float64, obj_val_L::Float64, u_vals::Vector{Float64}
    else
        inst.nTotalViewsF       = inst.nTotalViewsF ./ inst.Nsc
        inst.nOrganicViewsF     = inst.nOrganicViewsF ./ inst.Nsc
        inst.nForwardsF         = inst.nForwardsF ./ inst.Nsc
        inst.nTotalViewsL       = inst.nTotalViewsL ./ inst.Nsc
        inst.nOrganicViewsL     = inst.nOrganicViewsL ./ inst.Nsc
        inst.nForwardsL         = inst.nForwardsL ./ inst.Nsc
        return obj_val_F::Float64, obj_val_L::Float64, u_vals::Vector{Float64}
    end
end


# EVALUATE SEED SETS Ex-Post (on large set of scenarios Ω')
function evaluateObjectiveFunctionExPost( inst::instance, b::bfsdata, SL::Array{Int,1}, SF::Array{Int,1} )
    obj_val_F::Float64 = 0.0
    obj_val_L::Float64 = 0.0
    uF_resistant::Float64 = 0.0
    uL_resistant::Float64 = 0.0
    uF_organic::Float64 = 0.0
    uL_organic::Float64 = 0.0
    uF_total::Float64 = 0.0
    uL_total::Float64 = 0.0
    uF_forwards::Float64 = 0.0
    uL_forwards::Float64 = 0.0
    u_vals = zeros(Float64, inst.Nsc)
    inst.lgo_tmp  = [Tnodes[] for i=1:inst.nnodes]
    inst.lgfo_tmp = [Tnodes[] for i=1:inst.nnodes]
    inst.nTotalViewsF       = zeros(Float64, inst.nnodes)
    inst.nOrganicViewsF     = zeros(Float64, inst.nnodes)
    inst.nForwardsF         = zeros(Float64, inst.nnodes)
    inst.nTotalViewsL       = zeros(Float64, inst.nnodes)
    inst.nOrganicViewsL     = zeros(Float64, inst.nnodes)
    inst.nForwardsL         = zeros(Float64, inst.nnodes)

    for _ in 1:inst.Nsc
        clear!(b)
        if inst.kL == 0
            computeViewingValuesEvalIMP!( inst, b, SL, SF )
        else
            computeViewingValuesEvalCIMP!( inst, b, SL, SF )
        end


        for i in 1:inst.nnodes
            inst.nTotalViewsF[i]       += b.yv[i]
            inst.nOrganicViewsF[i]     += min(b.yv[i], 1)
            inst.nForwardsF[i]         += b.yf[i]
            inst.nTotalViewsL[i]       += b.xv[i]
            inst.nOrganicViewsL[i]     += min(b.xv[i], 1)
            inst.nForwardsL[i]         += b.xf[i]

            uF_resistant  += b.yv[i] / (b.xv[i] + b.yv[i] + inst.r[i])          #CAUTION: yv is also used for IMP here, thus xv=0
            uL_resistant  += (inst.leaderUtility * b.xv[i]) / ((inst.leaderUtility * b.xv[i]) + b.yv[i] + inst.r[i])          #CAUTION: yv is also used for IMP here, thus this is zero in the IMP
            uF_total      += b.yv[i]                                            #CAUTION: yv is also used for IMP here, thus xv=0
            uL_total      += b.xv[i]                                            #CAUTION: yv is also used for IMP here, thus this is zero in the IMP
            uF_forwards   += b.yf[i]                                            #CAUTION: yv is also used for IMP here, thus xv=0
            uL_forwards   += b.xf[i]
            #adjust viewings for organic case
            b.xv[i] > 1.0 ? b.xv[i] = 1.0 : nothing
            b.yv[i] > 1.0 ? b.yv[i] = 1.0 : nothing
            uF_organic      += b.yv[i]                                          #CAUTION: yv is also used for IMP here, thus xv=0
            uL_organic      += b.xv[i]                                          #CAUTION: yv is also used for IMP here, thus this is zero in the IMP
        end
    end

    uF_resistant    = uF_resistant * (1/inst.Nsc)
    uL_resistant    = uL_resistant * (1/inst.Nsc)
    uF_total        = uF_total * (1/inst.Nsc)
    uL_total        = uL_total * (1/inst.Nsc)
    uF_organic      = uF_organic * (1/inst.Nsc)
    uL_organic      = uL_organic * (1/inst.Nsc)
    uF_forwards     = uF_forwards * (1/inst.Nsc)
    uL_forwards     = uL_forwards * (1/inst.Nsc)
    inst.nTotalViewsF       = inst.nTotalViewsF ./ inst.Nsc
    inst.nOrganicViewsF     = inst.nOrganicViewsF ./ inst.Nsc
    inst.nForwardsF         = inst.nForwardsF ./ inst.Nsc
    inst.nTotalViewsL       = inst.nTotalViewsL ./ inst.Nsc
    inst.nOrganicViewsL     = inst.nOrganicViewsL ./ inst.Nsc
    inst.nForwardsL         = inst.nForwardsL ./ inst.Nsc


    return uF_resistant, uL_resistant, uF_total, uL_total, uF_organic, uL_organic, uF_forwards, uL_forwards
end



function computePropagationLengthIMP!( inst::instance, b::bfsdata, SL::Array{Int,1}, SF::Array{Int,1}, metric::String )
    avg_min  = 99999999.0
    avg_mean = 0.0
    avg_max  = 0.0
    length_count = 0

    # set seed sets to current level
    for i ∈ SF
        b.bvisited[i] = true
        push!(b.cur_level_F, i)
    end

    #create live graphs
    while !( isempty(b.cur_level_F) )
        #propagate from F
        @inbounds for i in b.cur_level_F
                        for t in inst.gout[i]
                            j, pf, pv = t[1], t[2], t[3]
                            if( pf == pv )                              #empirical pf > params["probV"]
                                if( rand( Bernoulli( pf ) ) == 1 )
                                    push!( inst.lgo_tmp[i], j )
                                    push!( inst.lgfo_tmp[i], j )
                                    if !b.bvisited[j]
                                        b.bvisited[j] = true
                                        push!(b.next_level_F, j)
                                    end
                                end
                            else
                                r = rand(Tprob)
                                if(r <= pf )
                                    push!( inst.lgo_tmp[i], j )
                                    push!( inst.lgfo_tmp[i], j )
                                    if !b.bvisited[j]
                                        b.bvisited[j] = true
                                        push!(b.next_level_F, j)
                                    end
                                elseif( pf < r <= pv )
                                    push!( inst.lgo_tmp[i], j )
                                else
                                    # do nothing
                                end
                            end
                        end
                    end
                empty!(b.cur_level_F)
                b.cur_level_F, b.next_level_F = b.next_level_F, b.cur_level_F
    end

    #compute propagation length
    myLiveGraph = inst.lgo_tmp
    if metric == "F" #forward max
        myLiveGraph = inst.lgfo_tmp
    end


    for s in SF
        fill!(b.bvisited, zero(Bool))
        nnodes = 0

        #get node ids /push seed node
        b.bvisited[s] = true
        push!(b.cur_level_F, s)

        length_count = 0
        while !( isempty(b.cur_level_F) )
            #propagate from F
            @inbounds for v in b.cur_level_F
                @inbounds for k in myLiveGraph[v]
                            if !b.bvisited[k]
                                b.bvisited[k] = true
                                push!(b.next_level_F, k)
                            end
                          end
                       end
            if !isempty(b.next_level_F)
                length_count += 1
            end
            empty!(b.cur_level_F)
            b.cur_level_F, b.next_level_F = b.next_level_F, b.cur_level_F
        end
#        println("l: $length_count")
        avg_mean += length_count
        avg_max < length_count ? avg_max = length_count : nothing
        avg_min > length_count ? avg_min = length_count : nothing

    end
    avg_mean = avg_mean / length(SF)


    for i=1:inst.nnodes
        empty!(inst.lgo_tmp[i])
        empty!(inst.lgfo_tmp[i])
    end

    return avg_min, avg_mean, avg_max
end


function evaluatePropagationLength( inst::instance, b::bfsdata, SL::Array{Int,1}, SF::Array{Int,1}, metric::String, wherefrom::String="OTHER_METHOD" )
    wi::Float64 = 0.0
    tempsum::Float64 = 0.0
    u_vals = zeros(Float64, inst.Nsc)
    inst.lgo_tmp  = [Tnodes[] for i=1:inst.nnodes]
    inst.lgfo_tmp = [Tnodes[] for i=1:inst.nnodes]
    avg_min  = 0.0
    avg_mean = 0.0
    avg_max  = 0.0
    tempsum_min  = 0.0
    tempsum_mean = 0.0
    tempsum_max  = 0.0

    for w = 1:inst.Nsc
        clear!(b)
        if inst.kL == 0
            avg_min, avg_mean, avg_max = computePropagationLengthIMP!( inst, b, SL, SF, metric )
            tempsum_min += avg_min
            tempsum_mean += avg_mean
            tempsum_max += avg_max
        else
            #CIMP case not implemented
        end
    end
    avg_min  = (1/inst.Nsc) * tempsum_min
    avg_mean = (1/inst.Nsc) * tempsum_mean
    avg_max  = (1/inst.Nsc) * tempsum_max

    return avg_min, avg_mean, avg_max
end
