#************************************************
#******** GET LEADERS SEED SET (IMP HEURISTIC)
#************************************************
# COMPUTE LEADERS SEED SET BASED ON LOCAL INDEX RANK (LIR) (cf., Liu2017)
function getSeedSetLeaderLIR( inst::instance )
    LI = Tuple[]                                                                #(node_id, LIR, outdegree_i)
    sizehint!(LI, inst.nnodes)
    S = Int[]                                                                   #leaders seed set
    sizehint!(S, inst.kL)

    for i=1:inst.nnodes
        li_val = 0
        Q = 0
        di = outdegree(inst.g, i)
        for j in outneighbors(inst.g, i)
            if ( (outdegree(inst.g, j) - di) > 0)
                li_val += 1
            else
                li_val += 0
            end
        end
        push!(LI, (i, li_val, di))
    end
    sort!(LI, by = x->(x[2],-x[3]))

    for k=1:inst.kL
        push!(S, LI[k][1] )
    end
    sort!(S)
    return S::Array{Int64,1}
end


#************************************************
#******** IMP HEURISTIC
#************************************************

# MARGINAL GAIN. NOTE: xv is abused as xf here for mtype == 1 (forwardings)
function heuristicMarginalGainIMP_all_models( inst::instance, b::bfsdata )
    S::Vector{Int} = Int[]                                                      #final set
    ReachedNodes::Vector{Vector{Tnodes}} = [Tnodes[] for w=1:inst.nsc]          #set of nodes that are reached in scenario w
    Rout::Matrix{Vector{Tnodes}} = [Tnodes[] for w=1:inst.nsc, j=1:inst.nnodes]
    cur_xv::Matrix{Float64} = zeros(Float64, inst.nnodes, inst.nsc)             #current value of x^vω[i]
    sizehint!(S, inst.kF)
    pq = PriorityQueue{Int,Float64}( Base.Order.Reverse )
    cur_best_node::Int = 0
    former_best_node::Int = -1
    cur_best_val::Float64 = -1
    tmp_sum::Float64 = 0.0
    neg_marginal_gain::Bool = false

    #compute marginal gain of each node (and sort in queue)
    for j=1:inst.nnodes
        tmp_sum = 0.0
        for w=1:inst.nsc
            clear!(b)
            for k ∈ inst.Rout[j,w]
                push!(Rout[w,j], k)
                if inst.mtype in [2,3,5]
                    for i ∈ inst.lgo[k,w]
                        b.xv[i]     += 1.0
                        cur_xv[i,w] += 1.0
                        b.bvisited[i] = true
                    end
                end
                if inst.mtype == 1
                    b.xv[k]     = 1.0
                    cur_xv[k,w] = 1.0
                    b.bvisited[k] = true
                end
            end
            if inst.mtype in [2,3,5]
                b.xv[j]   = 0.0
                cur_xv[j,w] = 0.0
                if inst.mtype == 5                  #patronage
                    for k ∈ 1:inst.nnodes
                        if b.bvisited[k]
                            tmp_sum += b.xv[k] / (b.xv[k] + inst.r[k])
                        end
                    end
                end
                if inst.mtype == 2                  #total views
                    for k ∈ 1:inst.nnodes
                        if b.bvisited[k]
                            tmp_sum += b.xv[k]
                        end
                    end
                end
                if inst.mtype == 3                  #organic reach
                    for k ∈ 1:inst.nnodes
                        b.xv[k] > 1.0 ? b.xv[k] = 1.0 : nothing
                        cur_xv[k,w] > 1.0 ? cur_xv[k,w] == 1.0 : nothing
                        if b.bvisited[k]
                            tmp_sum += b.xv[k]
                        end
                    end
                end
            end
            if inst.mtype == 1                      #forwarders
                for k ∈ 1:inst.nnodes
                    if b.bvisited[k]
                        tmp_sum += b.xv[k]
                    end
                end
            end
        end
        enqueue!(pq, j, tmp_sum )
    end

    #extract best node and add it to set
    cur_best_val = peek(pq)[2]
    cur_best_node = dequeue!(pq)
    push!(S, cur_best_node )
    for w=1:inst.nsc
        ReachedNodes[w] = Rout[w, cur_best_node]
    end


    #find node with best marginal gain (walk through queue until first element does not change)
    while( length(S) < inst.kF )
        while( cur_best_node != former_best_node )
            cur_best_node = peek(pq)[1]
            former_best_node = cur_best_node
            tmp_sum = 0.0
            neg_marginal_gain = false

            for w=1:inst.nsc
                clear!(b)
                #compute marginal viewings
                setdiff!( Rout[w, cur_best_node], ReachedNodes[w] )             #R^+_ω(cur_best_node) = R^+_ω(cur_best_node) ∖ R^+_ω(S)
                if inst.mtype in [2,3,5]
                    for k ∈ Rout[w, cur_best_node]
                        for i ∈ inst.lgo[k,w]
                            b.xv[i] += 1.0
                        end
                    end
                end
                if inst.mtype == 1
                    for k ∈ Rout[w, cur_best_node]
                        b.xv[k] = 1.0
                    end
                end

                if inst.mtype in [2,3,5]
                    b.xv[cur_best_node] = 0.0
                    for k ∈ S
                        b.xv[k] = 0.0
                    end
                end
                #compute marginal gain
                if inst.mtype == 5
                    for k ∈ 1:inst.nnodes
                        tmp_sum += ((cur_xv[k,w] + b.xv[k]) / (cur_xv[k,w] +b.xv[k] + inst.r[k])) - cur_best_val
                    end
                end
                if inst.mtype == 2                          #total
                    for k ∈ 1:inst.nnodes
                        tmp_sum += cur_xv[k,w] + b.xv[k] - cur_best_val
                    end
                end
                if inst.mtype == 3                          #organic
                    for k ∈ 1:inst.nnodes
                        b.xv[k] > 1.0 ? b.xv[k] = 1.0 : nothing
                        tmp_sum += cur_xv[k,w] + b.xv[k] - cur_best_val
                    end
                end
                if inst.mtype == 1                          #forwards
                    for k ∈ 1:inst.nnodes
                        tmp_sum += cur_xv[k,w] + b.xv[k] - cur_best_val
                    end
                end
            end
            pq[cur_best_node] = tmp_sum
            cur_best_node = peek(pq)[1]
            if tmp_sum < 0
                neg_marginal_gain = true
            end

            #check mem
            mb = get_mem_use()
            inst.memok = memOK(mb, inst.memlimit)
            if !inst.memok
                return Int[]                            #return an empty seed set since we could not solve MAR
            end
        end #end while cur_best_node != former_best_node

        push!( S, cur_best_node )

        #update viewings and current value
        tmp_sum = 0.0
        for w=1:inst.nsc
            #compute current viewings
            if inst.mtype in [2,3,5]
                for k ∈ Rout[w, cur_best_node]
                    for i ∈ inst.lgo[k,w]
                        cur_xv[i,w] += 1.0
                    end
                end
                for k ∈ S
                    cur_xv[k,w] = 0.0
                end
                if inst.mtype == 5
                    for k ∈ 1:inst.nnodes
                        tmp_sum += cur_xv[k,w] / (cur_xv[k,w] + inst.r[k])
                    end
                end
                if inst.mtype == 2
                    for k ∈ 1:inst.nnodes
                        tmp_sum += cur_xv[k,w]
                    end
                end
                if inst.mtype == 3
                    for k ∈ 1:inst.nnodes
                        cur_xv[k,w] > 1.0 ? cur_xv[k,w] == 1.0 : nothing
                        tmp_sum += cur_xv[k,w]
                    end
                end
            end
            if inst.mtype == 1
                for k ∈ Rout[w, cur_best_node]
                    cur_xv[k,w] += 1.0
                end
                for k ∈ 1:inst.nnodes
                    tmp_sum += cur_xv[k,w]
                end
            end
        end
        cur_best_val = tmp_sum

        #update current reached nodes
        for w=1:inst.nsc
            union!( ReachedNodes[w], Rout[w, cur_best_node] )
        end

        #dequeue
        dequeue_pair!(pq)
        cur_best_node = peek(pq)[1]
    end #end while |S| < k
    sort!(S)
    clear!(b)
    return S::Vector{Int}
end


# EXPECTED OUTDEGREE
function expected_outdegree_heuristic!( inst::instance, res::results )
    println("compute rankings...")
    inst.DEGRC  = Tuple{Tnodes,Float64}[]
    sizehint!(inst.DEGRC, inst.nnodes)

    #push tuples in vector
    for i in 1:inst.nnodes
        tempsum = 0.0
        for t in inst.gout[i]
            j, pf, pv = t[1], t[2], t[3]
            tempsum += pv
        end
        push!( inst.DEGRC, (i, tempsum))
    end

    #sort
    sort!(inst.DEGRC, by = x->(-x[2]) )

    #find ranks
    S = Int[]
    k=1
    ind=1
    while( k <= inst.kF )
        if( inst.DEGRC[ind][1] ∉ inst.SL )
            push!(S, inst.DEGRC[ind][1] )
            k += 1
            ind += 1
        else
            ind += 1
        end
    end
    sort!(S)

    println("e-outdegree SF:\t ", S)
    inst.SFDEG = S

    return
end


# TUNK RANK cf., https://thenoisychannel.com/2009/01/13/a-twitter-analog-to-pagerank
function computeTunkRank( inst::instance, itype::Int )
	tr = zeros(Float64, inst.nnodes)	#tunkranks
	ϵ  = 0.000001						#stopping criteria

	#initial tunkranks
	for i in 1:inst.nnodes
		tempsum = 0.0
		for tuple in inst.gout[i]
            j, pf, pv = tuple
            if itype == 3               #twitter
                if inst.nfriends[j] > 0
    				tempsum += 1 / inst.nfriends[j]
    			end
            else
                deg = indegree(inst.g, j)
                if deg > 0
    				tempsum += 1 / deg
    			end
            end
		end
		tr[i] = tempsum
	end

    #compute tunk ranks recursively with 100 iterations
	err = 0.0
	for _ in 1:100
		err = 0.0
		for i in 1:inst.nnodes
			last_tr_val = tr[i]
            cur_val_i   = 0.0
			for tuple in inst.gout[i]
                j, pf, pv = tuple
                if itype == 3       #twitter
                    if inst.nfriends[j] > 0
    					cur_val_i += (1 + pf * tr[j]) / inst.nfriends[j]
    				end
                else
                    deg = indegree(inst.g, j)
                    if deg > 0
    					cur_val_i += (1 + pf * tr[j]) / deg
    				end
    			end
            end
            tr[i] = cur_val_i
			err += abs( last_tr_val - tr[i] )
		end
		if err < inst.nnodes * ϵ
#            println("TunkRank converged after $_ iterations")
			return tr
		end
	end
	@warn "TuncRank did not converge with error $err"

	return tr
end

#REVERSE PAGE RANK (d=damping factor [0,1], nItr= number of recourse iterations)
function computeReversePageRank( inst::instance, d::Float64, nItr::Int, itype::Int )
	pr = zeros(Float64, inst.nnodes)	#pagerank
	ϵ  = 0.000001						#stopping criteria

	#initial tunkranks
	for i in 1:inst.nnodes
		tempsum = 0.0
		for tuple in inst.gout[i]
            j, pf, pv = tuple
			tempsum += (1-d) / inst.nnodes
		end
		pr[i] = tempsum
	end

    #compute reverse page ranks recursively with 100 iterations
	err = 0.0
	for _ in 1:nItr
		err = 0.0
		for i in 1:inst.nnodes
			last_tr_val = pr[i]
            cur_val_i   = 0.0
			for tuple in inst.gout[i]
                j, pf, pv = tuple
                if itype == 3 #twitter
                    if inst.nfriends[j] > 0
                        cur_val_i += (1-d) / inst.nnodes + d * pr[j] / inst.nfriends[j]
                    end
                else
                    indeg = indegree(inst.g, j)
                    if indeg > 0
                        cur_val_i += (1-d) / inst.nnodes + d * pr[j] / indeg
                    end
                end
			end
            pr[i] = cur_val_i
			err += abs( last_tr_val - pr[i] )
		end
		if err < inst.nnodes * ϵ
#            println("ReversePageRank converged after $_ iterations")
			return pr
		end
	end
	@warn "Reverse PageRank did not converge with error $err"

	return pr
end


"""
 OTHER HEURISTICS for (C)IMP 
 Node ranks are computed based on some criteria (e.g., betweeness, closeness, ),
 and the first k best nodes are chosen as seed node.
 Besides the computed ranks are compared with seed sets obtained with GBD or MAR, e.g.
 abs. ranks betweenness:  [(481, 3)] means that seed node 481 obtaind with GBD
 has a rank of 3 w.r.t. betweenness centrality.
"""
function computeRankings!( inst::instance, res::results, SF::Array{Int,1}, wherefrom::String, params::Dict, to::TimerOutput )
    println("compute rankings...")

    #store node id and correponding ranking criteria
    betweenness_ranks = Tuple{Int,Float64}[]
    degree_ranks      = Tuple{Int,Float64}[]
    closeness_ranks   = Tuple{Int,Float64}[]
    eoutdegree_ranks  = Tuple{Int,Float64}[]
    pagerank_ranks    = Tuple{Int,Float64}[]
    retweet_ranks     = Tuple{Int,Float64}[]        #retweets and mentions (cf., Leavitt2010)
    tunkrank_ranks    = Tuple{Int,Float64}[]
    sizehint!(betweenness_ranks, inst.nnodes)
    sizehint!(degree_ranks, inst.nnodes)
    sizehint!(closeness_ranks, inst.nnodes)
    sizehint!(eoutdegree_ranks, inst.nnodes)
    sizehint!(pagerank_ranks, inst.nnodes)
    sizehint!(retweet_ranks, inst.nnodes)
    sizehint!(tunkrank_ranks, inst.nnodes)

    #compute centrality measures (returns a vector with a value for each node)
    @timeit to "Hbc" bc = betweenness_centrality( inst.g )
    @timeit to "Hdc" dc = outdegree_centrality( inst.g )
    @timeit to "Hcc" cc = closeness_centrality( inst.g )

    #ReversePageRank
    @timeit to "Hpr" begin
        pr = computeReversePageRank( inst, 0.85, 100, params["itype"] )
        for i in 1:inst.nnodes
            push!(pagerank_ranks, (i, pr[i]) )
        end
        sort!(pagerank_ranks, by = x->(-x[2]) )
    end

    #TunkRank
    @timeit to "Htr" begin
        tr = computeTunkRank( inst, params["itype"] )
        for i in 1:inst.nnodes
            push!(tunkrank_ranks, (i, tr[i]) )
        end
        sort!(tunkrank_ranks, by = x->(-x[2]) )
    end

    #compute expected outdegrees
    @timeit to "Hed" begin
        for i in 1:inst.nnodes
            tempsum = 0.0
            for t in inst.gout[i]
                j, pf, pv = t[1], t[2], t[3]
                tempsum += pv
            end
            push!( eoutdegree_ranks, (i, tempsum))
        end
    end

    #push tuples in vector
    for i in 1:inst.nnodes
        push!( betweenness_ranks, (i, bc[i]) )
        push!( degree_ranks, (i, dc[i]) )
        if cc[i] == 0.0
            push!( closeness_ranks, (i, inst.nnodes) )         
        else
            push!( closeness_ranks, (i, cc[i]) )
        end
    end

    #sort by value
    @timeit to "Hbc" sort!(betweenness_ranks, by = x->(-x[2]))
    @timeit to "Hdc" sort!(degree_ranks, by = x->(-x[2]))
    @timeit to "Hcc" sort!(closeness_ranks, by = x->(-x[2]) )
    @timeit to "Hed" sort!(eoutdegree_ranks, by = x->(-x[2]) )


    #compute twitter ranks according to Leavitt2010
    if params["itype"] == 3
        #Leavitt2010
        @timeit to "Hrm" begin
            total_nretweets = sum( inst.nretweetsUserGets )
            total_nreplies  = sum( inst.nrepliesUserGets )
            total_nmentions = sum( inst.nmentionsUserGets )
            denom = total_nretweets + total_nreplies + total_nmentions
            for i in 1:inst.nnodes
                myrank = (inst.nretweetsUserGets[i] + inst.nrepliesUserGets[i] + inst.nmentionsUserGets[i]) / denom
                push!(retweet_ranks, (i, myrank) )
            end
            sort!(retweet_ranks, by = x->(-x[2]) )
        end
    end

    #find absolute position of seed nodes in rankings (k=seednode, i=rank), Vector has seedset size.
    function getAbsoluteRanksOfSeedNodes( myRanking )
        mySeedSetRanks = Tuple{Int,Int}[]
        for i in 1:inst.nnodes
            k, val = myRanking[i]
            if k ∈ SF
                push!(mySeedSetRanks, (k,i) )
            end
        end
        return mySeedSetRanks
    end

    #compute relative ranking (k=seednode, relRank=relative rank)
    function getRelativeRanksOfSeedNodes( myAbsRanks )
        mySeedSetRanksRel = Tuple{Int,Float64}[]
        for i in eachindex(myAbsRanks)
            k, absoluteRank = myAbsRanks[i]
            relativeRank    = absoluteRank / inst.nnodes
            push!(mySeedSetRanksRel, (k, relativeRank))
        end
        return mySeedSetRanksRel
    end

    #compute rankings
    if wherefrom == "BEN"
        res.abs_rank_bc = getAbsoluteRanksOfSeedNodes( betweenness_ranks )
        res.abs_rank_dc = getAbsoluteRanksOfSeedNodes( degree_ranks )
        res.abs_rank_cc = getAbsoluteRanksOfSeedNodes( closeness_ranks )
        res.abs_rank_ed = getAbsoluteRanksOfSeedNodes( eoutdegree_ranks )
        res.abs_rank_pr = getAbsoluteRanksOfSeedNodes( pagerank_ranks )
        res.abs_rank_tr = getAbsoluteRanksOfSeedNodes( tunkrank_ranks )
        res.rel_rank_bc = getRelativeRanksOfSeedNodes( res.abs_rank_bc )
        res.rel_rank_dc = getRelativeRanksOfSeedNodes( res.abs_rank_dc )
        res.rel_rank_cc = getRelativeRanksOfSeedNodes( res.abs_rank_cc)
        res.rel_rank_ed = getRelativeRanksOfSeedNodes( res.abs_rank_ed )
        res.rel_rank_pr = getRelativeRanksOfSeedNodes( res.abs_rank_pr)
        res.rel_rank_tr = getRelativeRanksOfSeedNodes( res.abs_rank_tr )
        if params["itype"] == 3
            res.abs_rank_rm = getAbsoluteRanksOfSeedNodes( retweet_ranks )
            res.rel_rank_rm = getRelativeRanksOfSeedNodes( res.abs_rank_rm )
        end

        if params["debug"] >= 3
            println("BEN:")
            println("abs. ranks betweenness:\t ",   res.abs_rank_bc)
            println("abs. ranks degree:\t\t ",      res.abs_rank_dc)
            println("abs. ranks closeness:\t ",     res.abs_rank_cc)
            println("abs. ranks e-outdegree:\t ",   res.abs_rank_ed)
            println("abs. ranks pagerank:\t ",      res.abs_rank_pr)
            println("abs. ranks retweet:\t\t ",     res.abs_rank_rm)
            println("abs. ranks tunkrank:\t ",      res.abs_rank_tr)
            println("rel. ranks betweenness:\t ",   res.rel_rank_bc)
            println("rel. ranks degree:\t\t ",      res.rel_rank_dc)
            println("rel. ranks closeness:\t ",     res.rel_rank_cc)
            println("rel. ranks e-outdegree:\t ",   res.rel_rank_ed)
            println("rel. ranks pagerank:\t ",      res.rel_rank_pr)
            println("rel. ranks retweet:\t\t ",     res.rel_rank_rm)
            println("rel. ranks tunkranks:\t ",     res.rel_rank_tr)
        end
    end
    if wherefrom == "MAR"
        res.abs_rank_bcMAR = getAbsoluteRanksOfSeedNodes( betweenness_ranks )
        res.abs_rank_dcMAR = getAbsoluteRanksOfSeedNodes( degree_ranks )
        res.abs_rank_ccMAR = getAbsoluteRanksOfSeedNodes( closeness_ranks )
        res.abs_rank_edMAR = getAbsoluteRanksOfSeedNodes( eoutdegree_ranks )
        res.abs_rank_prMAR = getAbsoluteRanksOfSeedNodes( pagerank_ranks )
        res.abs_rank_trMAR = getAbsoluteRanksOfSeedNodes( tunkrank_ranks )
        res.rel_rank_bcMAR = getRelativeRanksOfSeedNodes( res.abs_rank_bcMAR )
        res.rel_rank_dcMAR = getRelativeRanksOfSeedNodes( res.abs_rank_dcMAR )
        res.rel_rank_ccMAR = getRelativeRanksOfSeedNodes( res.abs_rank_ccMAR)
        res.rel_rank_edMAR = getRelativeRanksOfSeedNodes( res.abs_rank_edMAR )
        res.rel_rank_prMAR = getRelativeRanksOfSeedNodes( res.abs_rank_prMAR)
        res.rel_rank_trMAR = getRelativeRanksOfSeedNodes( res.abs_rank_trMAR )
        if params["itype"] == 3
            res.abs_rank_rmMAR = getAbsoluteRanksOfSeedNodes( retweet_ranks )
            res.rel_rank_rmMAR = getRelativeRanksOfSeedNodes( res.abs_rank_rmMAR )
        end
        if params["debug"] >= 2
            println("MAR")
            println("abs. ranks betweenness:\t ",   res.abs_rank_bcMAR)
            println("abs. ranks degree:\t\t ",      res.abs_rank_dcMAR)
            println("abs. ranks closeness:\t ",     res.abs_rank_ccMAR)
            println("abs. ranks e-outdegree:\t ",   res.abs_rank_edMAR)
            println("abs. ranks pagerank:\t ",      res.abs_rank_prMAR)
            println("abs. ranks retweet:\t\t ",     res.abs_rank_rmMAR)
            println("abs. ranks tunkrank:\t ",      res.abs_rank_trMAR)
            println("rel. ranks betweenness:\t ",   res.rel_rank_bcMAR)
            println("rel. ranks degree:\t\t ",      res.rel_rank_dcMAR)
            println("rel. ranks closeness:\t ",     res.rel_rank_ccMAR)
            println("rel. ranks e-outdegree:\t ",   res.rel_rank_edMAR)
            println("rel. ranks pagerank:\t ",      res.rel_rank_prMAR)
            println("rel. ranks retweet:\t\t ",     res.rel_rank_rmMAR)
            println("rel. ranks tunkrank:\t ",      res.rel_rank_trMAR)
        end
        
    end

    #get seed sets based on ranks
    #select the first best ones
    function getSeedSetsBasedOnRanking( myRanking )
        k::Int   = 1
        idx::Int = 1
        S = Int[]
        while( k <= inst.kF )
            node, val = myRanking[idx]
            if( node ∉ inst.SL )
                push!(S, node)
                k += 1
                idx += 1
             else
                idx += 1
             end
        end
        sort!(S)
        return S
    end
    @timeit to "Hbc" inst.SFbc = getSeedSetsBasedOnRanking( betweenness_ranks )
    @timeit to "Hcc" inst.SFcc = getSeedSetsBasedOnRanking( closeness_ranks )
    @timeit to "Hdc" inst.SFdc = getSeedSetsBasedOnRanking( degree_ranks )
    @timeit to "Hed" inst.SFed = getSeedSetsBasedOnRanking( eoutdegree_ranks )
    @timeit to "Hpr" inst.SFpr = getSeedSetsBasedOnRanking( pagerank_ranks )
    @timeit to "Htr" inst.SFtr = getSeedSetsBasedOnRanking( tunkrank_ranks )
    if params["itype"] == 3
        @timeit to "Hrm" inst.SFrm = getSeedSetsBasedOnRanking( retweet_ranks )
    end

    #get runtimes
    res.rtHbc = round( TimerOutputs.time(to["computeRankings"]["Hbc"]) / 10^9 ; digits=2 )
    res.rtHcc = round( TimerOutputs.time(to["computeRankings"]["Hcc"]) / 10^9 ; digits=2 )
    res.rtHdc = round( TimerOutputs.time(to["computeRankings"]["Hdc"]) / 10^9 ; digits=2 )
    res.rtHed = round( TimerOutputs.time(to["computeRankings"]["Hed"]) / 10^9 ; digits=2 )
    res.rtHpr = round( TimerOutputs.time(to["computeRankings"]["Hpr"]) / 10^9 ; digits=2 )
    res.rtHtr = round( TimerOutputs.time(to["computeRankings"]["Htr"]) / 10^9 ; digits=2 )
    if params["itype"] == 3
        res.rtHrm = round( TimerOutputs.time(to["computeRankings"]["Hrm"]) / 10^9 ; digits=2 )
    end

    return
end


#************************************************
#******** CIMP HEURISTICS
#************************************************
#NOTE: i abuse yv as yf here for mtype == 1
function heuristicMarginalGainCIMP_all_models( inst::instance, b::bfsdata )
    S::Vector{Int} = Int[]                                                      #final set
    ReachedNodes::Vector{Vector{Tnodes}} = [Tnodes[] for w=1:inst.nsc]          #set of nodes that are reached in scenario w
    Rout::Matrix{Vector{Tnodes}} = [Tnodes[] for w=1:inst.nsc, j=1:inst.nnodes]
    cur_yv::Matrix{Float64} = zeros(Float64, inst.nnodes, inst.nsc)             #current value of y^vω[i]
    sizehint!(S, inst.kF)
    pq = PriorityQueue{Int,Float64}( Base.Order.Reverse )
    cur_best_node::Int = 0
    former_best_node::Int = -1
    cur_best_val::Float64 = -1
    tmp_sum::Float64 = 0.0
    neg_marginal_gain::Bool = false

    #compute marginal gain of each node (and sort in queue)
    for j=1:inst.nnodes
        if j ∉ inst.SL
            tmp_sum = 0.0
            for w=1:inst.nsc
                clear!(b)
                #compute yv
                if inst.mtype in [2,3,5]
                    for k ∈ inst.Rout[j,w]
                        push!(Rout[w,j], k)
                        b.yfb[k] = true
                        for i ∈ inst.lgo[k,w]
                            b.yv[i]     += 1.0
                            cur_yv[i,w] += 1.0
                            b.bvisited[i] = true
                        end
                    end
                    b.yv[j]     = 0.0
                    cur_yv[j,w] = 0.0
                end
                if inst.mtype == 1
                    for k ∈ inst.Rout[j,w]
                        push!(Rout[w,j], k)
                        b.yfb[k] = true
                        b.yv[k]     = 1.0
                        cur_yv[k,w] = 1.0
                        b.bvisited[k] = true
                    end
                end


                #compute marginal gain
                if inst.mtype == 5
                    for k ∈ 1:inst.nnodes
                        if b.bvisited[k]
                            tmp_sum += b.yv[k] / (b.yv[k] + inst.rw[k,w])
                        end
                    end
                end
                if inst.mtype == 2
                    for k ∈ 1:inst.nnodes
                        if b.bvisited[k]
                            tmp_sum += b.yv[k]
                        end
                    end
                end
                if inst.mtype == 3
                    for k ∈ 1:inst.nnodes
                        b.yv[k] > 1.0 ? b.yv[k] = 1.0 : nothing
                        cur_yv[k,w] > 1.0 ? cur_yv[k,w] = 1.0 : nothing
                        if b.bvisited[k]
                            tmp_sum += b.yv[k]
                        end
                    end
                end
                if inst.mtype == 1
                    for k ∈ 1:inst.nnodes
                        if b.bvisited[k]
                            tmp_sum += b.yv[k]
                        end
                    end
                end

            end
            enqueue!(pq, j, tmp_sum )
        end
    end

    #extract best node and add it to set
    cur_best_val = peek(pq)[2]
    cur_best_node = dequeue!(pq)
    push!( S, cur_best_node )
    for w=1:inst.nsc
        ReachedNodes[w] = Rout[w, cur_best_node]
    end


    #find node with best marginal gain (walk through queue until first element does not change)
    while( length(S) < inst.kF )
        while( cur_best_node != former_best_node )
            cur_best_node = peek(pq)[1]
            former_best_node = cur_best_node
            tmp_sum = 0.0
            neg_marginal_gain = false

            for w=1:inst.nsc
                clear!(b)
                #compute marginal viewings
                setdiff!( Rout[w, cur_best_node], ReachedNodes[w] )             #R^+_ω(cur_best_node) = R^+_ω(cur_best_node) ∖ R^+_ω(S)
                if inst.mtype in [2,3,5]
                    for k ∈ Rout[w, cur_best_node]
                        for i ∈ inst.lgo[k,w]
                            b.yv[i] += 1.0
     #                        b.bvisited[i] = true
                        end
                    end
                    b.yv[cur_best_node] = 0.0
                    for k ∈ S
                        b.yv[k] = 0.0
                    end
                    for i in inst.SL
                        b.yv[i] = 0.0
                    end

                    #compute marginal gain
                    if inst.mtype == 5
                        for k ∈ 1:inst.nnodes
                            tmp_sum += ((cur_yv[k,w] + b.yv[k]) / (cur_yv[k,w] + b.yv[k] + inst.rw[k,w])) - cur_best_val
                        end
                    end
                    if inst.mtype == 2
                        for k ∈ 1:inst.nnodes
                            tmp_sum += cur_yv[k,w] + b.yv[k] - cur_best_val
                        end
                    end
                    if inst.mtype == 3
                        for k ∈ 1:inst.nnodes
                            b.yv[k] > 1.0 ? b.yv[k] = 1.0 : nothing
                            tmp_sum += cur_yv[k,w] + b.yv[k] - cur_best_val
                        end
                    end
                end
                if inst.mtype == 1
                    for k ∈ Rout[w, cur_best_node]
                        b.yv[k] = 1.0
                    end
                    for k ∈ 1:inst.nnodes
                        tmp_sum += cur_yv[k,w] + b.yv[k] - cur_best_val
                    end
                end

            end
            pq[cur_best_node] = tmp_sum
            cur_best_node = peek(pq)[1]
            if tmp_sum < 0
                neg_marginal_gain = true
            end
            #check mem
            mb = get_mem_use()
            inst.memok = memOK(mb, inst.memlimit)
            if !inst.memok
                return Int[]                                #return empty set if MAR needs to much mem (=no solution found)
            end
        end #end while cur_best_node != former_best_node

        #add to set
        push!( S, cur_best_node )

        if inst.inSL[cur_best_node]
            @warn "adding a node j ∈ L to F in MAR!!!"
        end

        #update viewings and current value
        tmp_sum = 0.0
        for w=1:inst.nsc
            clear!(b)
            #compute current viewings
            if inst.mtype == [2,3,5]
                for k ∈ Rout[w, cur_best_node]
                    for i ∈ inst.lgo[k,w]
                        cur_yv[i,w] += 1.0
                    end
                end
                for k ∈ S
                    cur_yv[k,w] = 0.0
                end

                for i in inst.SL
                    b.yv[i] = 0
                end

                if inst.mtype == 5
                    for k ∈ 1:inst.nnodes
                        tmp_sum += cur_yv[k,w] / (cur_yv[k,w] + inst.rw[k,w])
                    end
                end
                if inst.mtype == 2
                    for k ∈ 1:inst.nnodes
                        tmp_sum += cur_yv[k,w]
                    end
                end
                if inst.mtype == 3
                    for k ∈ 1:inst.nnodes
                        cur_yv[k,w] > 1.0 ? cur_yv[k,w] = 1.0 : nothing
                        tmp_sum += cur_yv[k,w]
                    end
                end
            end
            if inst.mtype == 1
                for k ∈ Rout[w, cur_best_node]
                    cur_yv[k,w] = 1.0
                end
                for k ∈ 1:inst.nnodes
                    tmp_sum += cur_yv[k,w]
                end
            end
        end
        cur_best_val = tmp_sum

        #update current reached nodes
        for w=1:inst.nsc
            union!( ReachedNodes[w], Rout[w, cur_best_node] )
        end

        #dequeue
        dequeue_pair!(pq)
        cur_best_node = peek(pq)[1]
    end #while |S| < k
    sort!(S)
    clear!(b)
    return S::Vector{Int}
end