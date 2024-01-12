# INSTANCE READER
#include("var_declarations.jl")


# READ ARTIFICIAL INSTANCES
function read_Artificial_Instance!( inst::instance, params::Dict )
    print( "Reading instance..." )

    #set seed
    Random.seed!( params["rseed"] )

    #open file
    f = open( joinpath(inst.ipath, inst.ifile) )

    #read all lines as a string
    lines = readlines(f)

    # nnodes, nedges (written in line nr. 4)
    inst.nnodes = parse( Tnodes, split(lines[4])[1] )
    inst.narcs  = parse( Tarcs, split(lines[4])[2] )
    inst.gout = [Tuple{Tnodes,Tnodes}[] for i=1:inst.nnodes]

    # check for errors
    @assert( inst.narcs <= typemax( Tarcs ), "narcs > typemax(Tnodes)" )
    @assert( inst.nnodes <= typemax( Tnodes ), "nnodes > typemax(Tnodes)" )
    @assert( inst.nnodes > 0, "invalid nnodes" )
    @assert( inst.narcs > 0, "invalid narcs" )

    # create graph and read graph structure
    inst.g = SimpleDiGraph( inst.nnodes )
    sizehint!( inst.edge, inst.narcs )

    for l in 6:length( lines )
        i = parse( Tnodes, split(lines[l])[1] ) + 1            # add +1 since instance starts with node id 0, but we need to start with 1
        j = parse( Tnodes, split(lines[l])[2] ) + 1
        @assert( i != j, "graph contains self-loop at line $l" )
        @assert( i <= inst.nnodes || j <= inst.nnodes, "instance contains a vertex > nnodes at line $l")
        if (! add_edge!(inst.g, i, j) )
            error( "instance_reader.jl: arc already exists: ($i, $j)" )
        end
    end

    #delete singletons
    for i ∈ inst.nnodes:-1:1
        degree(inst.g, i) == 0 ? rem_vertex!( inst.g, i ) : nothing
    end
    inst.nnodes = nv( inst.g )

    #add my structs
    for a ∈ edges( inst.g )
        i, j = src(a), dst(a)
		pf = 0.01
        pv = max( pf, params["probV"] )
        @assert( pf <= pv, "pf > pv" )
        push!( inst.edge, (i, j , pf , pv) )
        push!( inst.pF, pf )
        push!( inst.gout[i], (j, pf, pv) )
    end

    # check for errors
    @assert( inst.narcs  == length( lines ) - 5, "narcs in instance file does not equal the number of read lines" )
    @assert( inst.narcs  == ne( inst.g ), "narcs in instance file does not equal ne(inst.g)" )
    @assert( inst.nnodes == nv( inst.g ), "nnodes in instance file does not equal nv(inst.g)" )

    #close file
    close(f)
    print("Done!\n")
    return nothing
end


# READ SNAP INSTANCES
function read_Snap_Instance!( inst::instance, params::Dict )
    Random.seed!( params["rseed"] )
    stripChar  = (s, r) -> replace(s, Regex("[$r]") => " ")


    #get forwarding probabilities extracted from twitter (used for boot strapping)
	filename = ""
    prob_filename = "test_probabilities.csv"
#    prob_filename = "tw-inststatsedge_rt_by_ot.csv"
    filename = joinpath(params["root_path"], "data", "edgestats", prob_filename)

    f                       = open( filename )
    lines::Vector{String}   = readlines(f)
    tw_probs::Vector{Tprob} = zeros(Tprob, length(lines)-1)
    sampleSize::Int         = length(tw_probs)
	for l in 2:length( lines )
		myline = stripChar( lines[l], ";")
		tw_probs[l-1] = parse( Tprob, lines[l] )
	end
    close(f)

    #open file
    f = open( joinpath(inst.ipath, inst.ifile) )

    #read all lines as a string
    lines = readlines(f)

    # nnodes, nedges (written in line nr. 3)
    inst.nnodes = parse( Tnodes, split(lines[3])[1] )
    inst.narcs  = parse( Tarcs, split(lines[3])[2] )
    nnodes_with_singletons = deepcopy( inst.nnodes )
    sizehint!( inst.edge, inst.narcs )

    # consitency check
    @assert( inst.narcs <= typemax( Tarcs ), "narcs > integer range" )
    @assert( inst.nnodes <= typemax( Tnodes ), "nnodes > integer range" )

    # create graph and read graph structure
    inst.g = SimpleDiGraph( inst.nnodes )


    for l in 5:length( lines )
#        println("line $l")
        i  = parse( Tnodes, split(lines[l])[1] ) + 1            # add +1 since instance starts with node id 0, but we need to start with 1
        j  = parse( Tnodes, split(lines[l])[2] ) + 1
        ne = parse( Int, split(lines[l])[3] )
        @assert( i <= inst.nnodes || j <= inst.nnodes, "instance contains a vertex > nnodes at line $l")

        #no self-edges
        if( i == j )
            inst.narcs -= 1
            continue
        end

        # no multi-arcs
        if( !has_edge(inst.g, i, j) && i != j )
            if(!add_edge!(inst.g, i, j))
                println("did not add edge $i $j at line $l")
                inst.narcs -=1
            else
                #nothing
            end
        else
            inst.narcs -=1
        end
    end

    #delete singletons
    for i ∈ inst.nnodes:-1:1
        degree(inst.g, i) == 0 ? rem_vertex!( inst.g, i ) : nothing
    end
    inst.nnodes = nv( inst.g )
    nremoved_vertices = nnodes_with_singletons - inst.nnodes
    println("Removed $nremoved_vertices singletons")


    # add edge data to instance struct
    inst.gout = [Tuple{Tnodes,Tnodes}[] for i=1:inst.nnodes]
#    inst.EiDeg = zeros(Tprob, inst.nnodes)
#    inst.EoDeg = zeros(Tprob, inst.nnodes)
#    params["computeInstanceStats"] ? inst.gin = [Tuple{Tnodes,Tnodes}[] for i=1:inst.nnodes] : nothing
    for a ∈ edges( inst.g )
        i, j = src(a), dst(a)
        ppos = rand(1:sampleSize)
        pf = tw_probs[ppos]
        pv = max( pf, params["probV"] )
        @assert( pf <= pv, "pf > pv" )
        push!( inst.edge, (i, j , pf , pv) )
        push!( inst.pF, pf )
        push!( inst.gout[i], (j, pf, pv) )
    end

    # consitency check
    na = inst.narcs
    nedges = ne( inst.g )
    @assert( inst.narcs  == ne( inst.g ), "narcs in instance file does not equal ne(inst.g) - $na != $nedges" )
    @assert( inst.nnodes == nv( inst.g ), "nnodes in instance file does not equal nv(inst.g)" )
	if !is_connected( inst.g )
		@warn "graph is not connected"
	end

    #close file
    close(f)
    print("Done!\n")
    println("nnodes=", inst.nnodes)
    println("narcs=", inst.narcs)
    return nothing
end


# READ TWITTER INSTANCES
function read_Twitter_Instance_MetaGraphs!( inst::instance, params::Dict )
    println("Reading Twitter instance...")
    #open file
    f = open( joinpath(inst.ipath, inst.ifile) )

    #read all lines as a string
    lines = readlines(f)

    # nnodes, nedges (written in line nr. 3)
    inst.nnodes = parse( Tnodes, split(lines[4])[1] )
    inst.narcs  = parse( Tarcs, split(lines[4])[2] )
#    inst.gout = [Tuple{Tnodes,Tnodes}[] for i=1:inst.nnodes]
#    params["computeInstanceStats"] ? inst.gin = [Tuple{Tnodes,Tnodes}[] for i=1:inst.nnodes] : nothing

    #sizehints
    sizehint!( inst.tname, inst.nnodes )
    sizehint!( inst.nfriends, inst.nnodes )
    sizehint!( inst.nfollowers, inst.nnodes )
    sizehint!( inst.notweets, inst.nnodes )
    sizehint!( inst.nretweets, inst.nnodes )
    sizehint!( inst.nreplies, inst.nnodes )
    sizehint!( inst.nlikes, inst.nnodes )
	sizehint!( inst.nretweetsUserGets, inst.nnodes )
    sizehint!( inst.nrepliesUserGets, inst.nnodes )
    sizehint!( inst.nmentionsUserGets, inst.nnodes )
    sizehint!( inst.edge, inst.narcs )


    # consitency check
    @assert( inst.narcs <= typemax( Tarcs ), "narcs > integer range" )
    @assert( inst.nnodes <= typemax( Tnodes ), "nnodes > integer range" )

	function check_str2(a) 	#if a twitter name as a blank, just use the first part
    	return tryparse(Int, a) !== nothing
	end

    # create metadigraph and read graph structure
    mg = MetaDiGraph( inst.nnodes, 0.0 )

    # read instance data from file
    reached_edge_list = false
    for l in 6:length( lines )
        if( split(lines[l])[1] == "#" )
            reached_edge_list = true
            continue
        end

        if( !reached_edge_list ) # read data for each twitter user profile in the instance
			set_prop!( mg, l-5, :tname, split(lines[l])[3] )        # user profile name)
			k=4
			while !check_str2( split(lines[l])[k] )
				if !check_str2( split(lines[l])[k] )
					k += 1
				end
			end
			set_prop!( mg, l-5, :nfriends, parse(Int, split(lines[l])[k]) )         # of friends
			set_prop!( mg, l-5, :nfollowers, parse(Int, split(lines[l])[k+1]) )     # of followers
			set_prop!( mg, l-5, :notweets, parse(Int, split(lines[l])[k+2]) )       # of original tweets
			set_prop!( mg, l-5, :nretweets, parse(Int, split(lines[l])[k+3]) )      # of retweets (i.e., user is retweeting)
			set_prop!( mg, l-5, :nreplies, parse(Int, split(lines[l])[k+4]) )       # of replies (i.e., user replies to someone)
			set_prop!( mg, l-5, :nlikes, parse(Int, split(lines[l])[k+5]) )         # of likes (i.e., user likes tweets of someone)
			set_prop!( mg, l-5, :nretweetsUserGets, 0 )                             # of retweets (i.e., user gets retweeted)
			set_prop!( mg, l-5, :nmentionsUserGets, 0 )                             # of mentions (i.e., user gets mentioned)
			set_prop!( mg, l-5, :nrepliesUserGets, 0 )                              # of replies (i.e., user gets replied)

        else #edge list
            addarc = false

			#parse data
            i = parse( Tnodes, split(lines[l])[2] ) + 1         #influencing user    # add +1 since instance starts with node id 0, but we need to start with 1
            j = parse( Tnodes, split(lines[l])[3] ) + 1         #influenced user
            nretweets_j = parse( Int, split(lines[l])[4] )      # number of retweets of j of influencer i
			nmentions_j = parse( Int, split(lines[l])[5] )      # number of answers of j to influencer i
            nanswers_j = parse( Int, split(lines[l])[6] )       # number of answers of j to influencer i

			#get current data of node i
			cur_retweets_i_gets = get_prop(mg, i, :nretweetsUserGets)
			cur_mentions_i_gets = get_prop(mg, i, :nmentionsUserGets)
			cur_replies_i_gets  = get_prop(mg, i, :nrepliesUserGets)


			notweets_i = get_prop(mg, i, :notweets)
			nretweets_i = get_prop(mg, i, :nretweets)
            tot_notweets_j  = get_prop(mg, j, :notweets)
			tot_nretweets_j = get_prop(mg, j, :nretweets)
			tot_nreplies_j = get_prop(mg, j, :nreplies)
            if( (tot_notweets_j + tot_nretweets_j + tot_nreplies_j) > 0 && (nretweets_j + nanswers_j > 0) ) # submitted version
                addarc = true
            else
                inst.narcs -=1
                continue
            end

            if( addarc )    # only add arc if influence is exerted
				# update influence
				set_prop!( mg, i, :nretweetsUserGets, cur_retweets_i_gets + nretweets_j)
				set_prop!( mg, i, :nmentionsUserGets, cur_mentions_i_gets + nmentions_j )
				set_prop!( mg, i, :nrepliesUserGets, cur_replies_i_gets + nanswers_j )
				tot_notweets_j  = get_prop(mg, j, :notweets)
				tot_nretweets_j = get_prop(mg, j, :nretweets)
				tot_nreplies_j = get_prop(mg, j, :nreplies)
				influence_probability::Tprob = min( (nretweets_j + nanswers_j) / (tot_notweets_j + tot_nretweets_j + tot_nreplies_j), 1.0) # submitted version

                # consitency check
                @assert( 0.0 <= influence_probability <= 1, "not 0 <= $influence_probability <=1 at line $l" )
                @assert( i != j, "graph contains self-loop at line $l" )
                @assert( i <= inst.nnodes || j <= inst.nnodes, "instance contains a vertex > nnodes at line $l")

                # no multi-arcs, or self-loops
                if( !has_edge(mg, i, j) && i != j )
                    if(!add_edge!(mg, i, j))
                        println("did not add edge $i $j at line $l")
                        inst.narcs -=1
                    else
						set_prop!(mg, i, j, :weight, influence_probability)
                    end
                else
                    inst.narcs -=1
                end
            end
        end
    end

	#delete singletons
	for i in nv(mg):-1:1
		if degree(mg,i) == 0
			rem_vertex!(mg, i)
		end
	end

	#build graph and set instance struct values
	inst.nnodes = nv(mg)
	inst.narcs  = ne(mg)
	inst.g = SimpleDiGraph( inst.nnodes )
	pv = params["probV"]
	inst.gout = [Tuple{Tnodes,Tnodes}[] for i=1:inst.nnodes]
	for i in 1:inst.nnodes
		push!( inst.tname, get_prop(mg, i, :tname) )
		push!( inst.nfriends, get_prop(mg, i, :nfriends) )
		push!( inst.nfollowers, get_prop(mg, i, :nfollowers) )
		push!( inst.notweets, get_prop(mg, i, :notweets) )
		push!( inst.nretweets, get_prop(mg, i, :nretweets) )
		push!( inst.nreplies, get_prop(mg, i, :nreplies) )
		push!( inst.nlikes, get_prop(mg, i, :nlikes) )
		push!( inst.nretweetsUserGets, get_prop(mg, i, :nretweetsUserGets) )
		push!( inst.nrepliesUserGets, get_prop(mg, i, :nrepliesUserGets) )
		push!( inst.nmentionsUserGets, get_prop(mg, i, :nmentionsUserGets) )
		for j in outneighbors(mg, i)
			add_edge!(inst.g, i, j)
			pf = get_prop(mg, i, j, :weight)
			pv = max( pf, pv )
			push!( inst.edge, (i, j , pf, pv) )
			push!( inst.pF, pf )
			push!( inst.gout[i], (j, pf, pv) )
		end
	end
	mg = nothing

    # consitency check
    inst.narcs = ne( inst.g )
    nedges = ne( inst.g )
    @assert( inst.narcs  == length(inst.edge), "narcs in instance file does not equal ne(inst.g) - $na != $nedges" )
    @assert( inst.nnodes == nv( inst.g ), "nnodes in instance file does not equal nv(inst.g)" )
	if !is_connected( inst.g )
		@warn "graph is not connected"
	end

    #close file
    close(f)
    print("Done!\n")
    println("nnodes=", inst.nnodes)
    println("narcs=", inst.narcs)

    return nothing
end


# READ LEADER SEED SET
function read_SSL!( inst::instance, params::Dict )
    print( "Reading SSL..." )
	mykL = "foo"
	if inst.kL < 10
		mykL = "0" * string(inst.kL)
	else
		mykL = string(inst.kL)
	end

    ssl_fname = "SSL-LOG_BEN_" * inst.ifile * "_k=" * mykL * "_m=5_r=0.5_w=100"     #this is the setting used in the paper
    fname = joinpath(params["root_path"], "data", "leader_seed_sets", ssl_fname)

    println("from: ", fname)

    if inst.kL == 0         #empty leader
        inst.SL = Int[]
    elseif isfile( fname )  #file exists
        f = open( fname )

        #read all lines as a string
        lines = readlines(f)

        # nnodes, nedges (written in line nr. 3)
        for l=1:length(lines)
            kL = parse( Int, split(lines[l])[1] )

            if( kL == inst.kL )
                for j = 2:inst.kL + 1
                    seed_node = parse( Int, split(lines[l])[j] )
                    push!(inst.SL, seed_node)
                end
            end
        end
        close(f)
    else
        @warn("SSL file does not exist (generating LIR seedset)")
        inst.SL = getSeedSetLeaderLIR( inst )
    end

    println("####SL:", inst.SL)
end
