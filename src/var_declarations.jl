# DECLARATION OF VARIABLES AND CONSTRUCTORS

#TYPE DEFINITONS
const Tnodes = UInt16
const Tarcs  = UInt32
const Tprob  = Float32


# INSTANCE STRUCT
mutable struct instance
    memok::Bool                         # true if memory is ok
    memlimit::Int                       # memory limit in mb
    ipath::String                       # instance path (e.g., /my/path/to/instance/)
    ifile::String                       # instance file (e.g., instance.txt)
    mtype::Int                          # model type
    g::SimpleDiGraph                    # instance graph G =( V, A )
    gout::Vector{Vector{Tuple{Tnodes,Tprob,Tprob}}}  # gout[i] -> tuples (j, pf, pv) outneighbors of node i in G^ω
    gin::Vector{Vector{Tuple{Tnodes,Tprob,Tprob}}}   # gout[i] -> inneighbors of node i in G
    edge::Vector{Tuple{Tnodes,Tnodes,Tprob,Tprob}}   # Tuple (i,j,pf,pv) (src,dst,prob forward, prob viewing)
    EoDeg::Vector{Tprob}                # expected out degree
    EiDeg::Vector{Tprob}                # expected out degree
    nnodes::Int                         # |V|
    narcs::Int                          # |A|
    nlivearcs::Vector{Int}              # number of active live arcs per scenario
    v::Array{UInt8, 1}                  # node weights
    r::Vector{Float64}                  # resistence values r_i
    rw::Matrix{Float64}                 # resistence values r^ω_i or rw[i,w]
    arcCounter::Array{Int,1}            # counts how often an arc was diced in live graphs
    pF::Array{Tprob,1}                  # probabilies for live-arc graphs |A|
    probF::Tprob                        # probability for live-arc graphs (if set constant)
    pvrange::Tuple{Tprob,Tprob}         # (pv_lower_bound, pv_upper_bound) viewing bounds
    l_LPbound::Float64                  # last LP relaxation value
    c_LPbound::Float64                  # current LP relaxation value
    l_expnode::Int                      # # of explored nodes in B&B tree
    c_expnode::Int                      # current eplored node
    kF::Int                             # seed-set size Follower
    kL::Int                             # seed-set size Leader
    leaderUtility::Int                  # factor for scaling leader utility
    SF::Array{Int,1}                    # seed-set nodes of Follower (from solution!)
    inSF::BitArray{1}                   # true if i ∈ SF, false otherwise (from solution!)
    SL::Array{Int,1}                    # seed-set nodes of Leader
    inSL::BitArray{1}                   # true if i ∈ SL, false otherwise
    SFSAA::Vector{Vector{Int}}          # seed sets obtained by sample average approximation
    SFMAR::Vector{Int}                  # seed set Marginal Gain heuristic
    SFMARSAA::Vector{Vector{Int}}       # seed sets obtained by sample average approximation
    SFDEG::Vector{Int}                  # seed set expected out-degree heuristic
    SFbc::Vector{Int}                   # seed set betweenness centralityy heuristic
    SFcc::Vector{Int}                   # seed set closeness centrality heuristic
    SFdc::Vector{Int}                   # seed set degree centrality heuristic
    SFed::Vector{Int}                   # seed set expected out-degree heuristic
    SFpr::Vector{Int}                   # seed set page rank
    SFrm::Vector{Int}                   # seed set retweets / mentions etc
    SFtr::Vector{Int}                   # seed set tunk rank
    SFmx::Vector{Int}                   # seed set from matrix powering heuristic
#    SFBIN::Vector{Int}                 # seed set Marginal Gain heuristic
    DEGRC::Vector{Tuple{Tnodes,Float64}}# compititive rank best indivuduals (node_id, rank=expected nr of reached nodes)
    nsc::Int                            # |Ω|
    Nsc::Int                            # |Ω'|, i.e., for evaluating the function
    nSSAItr::Int                        # number of sample average approximation iteratrion
    TLmin::Array{Int, 2}                # T_min(ω)_i  earliest activation time step Leader
    Rout::Matrix{Vector{Tnodes}}        # reachable set: Rout[i,w]: nodes that are activated by node i in scenario w
    Rin::Matrix{Vector{Tnodes}}         # reachableby set: Rin[i,w]: nodes that can activate by node i in scenario w
    RoutL::Vector{Vector{Tnodes}}       # Set of nodes the leader can activate if F=∅
    Singleton::BitArray{2}              # Singleton[i,w] = true indicates that node i is a singleton in scenario ω
    isSingleton::Array{Bool,1}          # isSingleton[i] = true if node i is a singleton in AT LEAST one scenario
    Star::BitArray{2}                   # Star[i,w] = true indicates that node i is a singleton in scenario ω
    MNC::Matrix{Tnodes}                 # F^ω(L,{j}), i.e., marginal node contributions ∀i ∈ V \L, ∀ω ∈ Ω (MNC[i,w])
    ireachesL::Array{Bool, 2}           # used for Bharati heuristic ireachesL[i,w] = true if i reaches a leader path in scenario w
    tmpContrib::Vector{Int}             # temporary values of contributions per scenario
    nActivations::Array{Float64, 2}     # Expected number of activations per node nActivations[by leader, by follower] ; length =|V|
    lgo::Matrix{Vector{Tnodes}}         # live-graph lgo[i,w] = (viewing + forwarding) Set of all outneighbors of node i in scenario w
    lgi::Matrix{Vector{Tnodes}}         # live-graph lgi[i,w] = (viewing + forwarding) Set of all inneighbors of node i in scenario w
    lgfo::Matrix{Vector{Tnodes}}        # live-graph lgo[i,w] = (only forwarding) Set of all outneighbors of node i in scenario w
    lgfi::Matrix{Vector{Tnodes}}        # live-graph lgi[i,w] = (only forwarding) Set of all inneighbors of node i in scenario w
    lgo_tmp::Vector{Vector{Tnodes}}     # used for evaluating obj function
    lgfo_tmp::Vector{Vector{Tnodes}}    # used for evaluating obj function
    indeg::Matrix{Tnodes}               # indegree indeg[i,w] = indegree in of node i in scneario ω
    tname::Array{String,1}              # twitter name
    nfriends::Array{Int64,1}            # twitter #friends
    nfollowers::Array{Int64,1}          # twitter #followers
    notweets::Array{Int64,1}            # twitter #otweets
    nretweets::Array{Int64,1}           # twitter #how often a user retweets something
    nreplies::Array{Int64,1}            # twitter #how often a user replies to someone
    nlikes::Array{Int64,1}              # twitter #likes a user receives
    nretweetsUserGets::Array{Int64,1}   # twitter #how often a user is retweeted
    nrepliesUserGets::Array{Int64,1}    # twitter #how often a user gets replies
    nmentionsUserGets::Array{Int64,1}   # twitter #how often a user is mentioned
    nTotalViewsF::Vector{Float64}       # expected number of total views Follower
    nOrganicViewsF::Vector{Float64}     # expected number of organic views Follower
    nForwardsF::Vector{Float64}         # expected number of tweets forwarded triggered by the Follower
    nTotalViewsL::Vector{Float64}       # expected number of total views Follower
    nOrganicViewsL::Vector{Float64}     # expected number of organic views Follower
    nForwardsL::Vector{Float64}         # expected number of tweets forwarded triggered by the Follower
    hashval::UInt64                     # hash value of [ifile, nsc, kL, kF]
end


# CREATE EMPTY INSTANCE
function createEmptyInstance( params::Dict )
    ipath::String = params["ipath"]
    ifile::String = params["ifile"]
    mtype::Int    = params["mtype"]
    kF::Int       = params["kF"]
    kL::Int       = params["kL"]
    nsc::Int      = params["nscenarios"]
    Nsc::Int      = params["Nscenarios"]
    nSSAItr::Int  = params["nSSAItr"]
    leaderUtility::Int  = params["leaderUtility"]

    return instance( true,
                     0,
                     ipath,
                     ifile,
                     mtype,
                     SimpleDiGraph(2),
                     Vector{Tuple{Tnodes,Tprob,Tprob}}[],
                     Vector{Tuple{Tnodes,Tprob,Tprob}}[],
                     Tuple{Tnodes,Tnodes,Tprob,Tprob}[],
                     Tprob[], Tprob[],
#                     "X",
                     -1,-1,
                     Int[],
                     UInt8[],
                     Float64[],
                     zeros(Float64,2,2),
                     Int[],
                     Float64[], 0.1, (0.3, 0.7),
                     typemax(Float64), 0.0, -1, -1,
                     kF, kL, leaderUtility,
                     Int[], falses(2),
                     Int[], falses(2),
#                     Int[], Int[], Int[], Int[],
                     Vector{Int}[],
                     Int[],
                     Vector{Int}[],
                     Int[], Int[], Int[], Int[], Int[], Int[], Int[], Int[], Int[], #seed sets F
                     Tuple{Tnodes,Float64}[],
#                     Tuple[],
#                     Tuple{Tnodes,Float64}[],
#                     -1, -1, -1, -1,
                     nsc, Nsc, nSSAItr,
#                     -1,
#                     zeros(Int, 2, 2),
#                     spzeros(Int, 2, 2),
                     zeros(Int, 2, 2),
#                     falses(2,2,2),
                     Matrix{Vector{Tnodes}}(undef,0,0),
                     Matrix{Vector{Tnodes}}(undef,0,0),
                     Vector{Tnodes}[],
                     falses(2,2),
                     zeros(Bool, 1),
                     falses(2,2),
                     zeros(Tnodes, 2, 2),
                     falses(2,2),
                     Int[],
                     zeros(Float64, 2, 2),
                     Matrix{Vector{Tnodes}}(undef,0,0),
                     Matrix{Vector{Tnodes}}(undef,0,0),
                     Matrix{Vector{Tnodes}}(undef,0,0),
                     Matrix{Vector{Tnodes}}(undef,0,0),
                     Vector{Tnodes}[],
                     Vector{Tnodes}[],
                     zeros(Tnodes, 2, 2),
#                     Float64[],
                     String[], Int[], Int[], Int[], Int[], Int[], Int[], Int[], Int[], Int[],
                     Float64[], Float64[], Float64[], Float64[], Float64[], Float64[],
                     0 )::instance
end


# RESULTS STRUCT
mutable struct results
    opath::String                       # output path (e.g., /my/path/to/results/)
    ofile::String                       # ouput filename (e.g., output.txt)
    cpxStatus::String                   # cplex status
    cpxStatusLP::String                 # cplex status (LP RELAXATION)
    bestboundSAA::Vector{Float64}       # store the best bounds found by cplex
    bestboundMARSAA::Vector{Float64}    # store the best bounds found by MAR heuristic
    bestbound::Float64                  # average upper bound for obj val
    obj_F::Float64                      # average objective value (model), i.e., Follower
    obj_L::Float64                      # objective value of Leader (computed in solution checker)
    obj_F_solcheck::Float64             # objective value obtained from solutionc checker
    obj_F_eval::Float64                 # objective value Follower after evaluation
    obj_L_eval::Float64                 # objective value of Leader after evaluation
    obj_FMAR::Float64                   # objective value marginal gain heuristic
    obj_FMAR_eval::Float64              # objective value marginal gain heuristic, evalueted
    obj_LMAR_eval::Float64              # objective value leader marginal gain heuristic, evaluated
    obj_FDEG::Float64                   # objective value follower expected outdegree heuristic
    obj_Fbc_eval::Float64               # objective value betweenness centrality eval
    obj_Lbc_eval::Float64               # objective value betweenness centrality eval
    obj_Fcc_eval::Float64               # objective value closeness centrality eval
    obj_Lcc_eval::Float64               # objective value closeness centrality eval
    obj_Fdc_eval::Float64               # objective value degree centrality eval
    obj_Ldc_eval::Float64               # objective value degree centrality eval
    obj_Fed_eval::Float64               # objective value expected out centrality eval
    obj_Led_eval::Float64               # objective value expected out centrality eval
    obj_Fpr_eval::Float64               # objective value pagerank
    obj_Lpr_eval::Float64               # objective value pagerank
    obj_Frm_eval::Float64               # objective value retweets / mentions
    obj_Lrm_eval::Float64               # objective value retweets / mentions
    obj_Ftr_eval::Float64               # objective value tunk rank eval
    obj_Ltr_eval::Float64               # objective value tunk rank eval
    obj_Fmx_eval::Float64               # objective value matrix powering heuristic
    LPrelax::Float64                    # LP relaxation
    LPgap::Float64                      # MIPgap in root relaxation [%]
    gap::Float64                        # optimiality gap
    gap_approx_abs::Float64             # estimated absolute approximation gap
    gap_approx_rel::Float64             # estimated relative approximation gap in percent!!
    gap_approx_relMAR::Float64          # estimated relative approximation gap in percent!!
    fcut_error::Float64                 # relative error between solution checker result and cplex
    lcdl::Array{Float64,1}              # lower confidence level of solution 1=90%, 2=95%, 3=99%
    ucdl::Array{Float64,1}              # upper confidence level of solution
    nBENCUTS::Int                       # number of Benders Cuts (calls of function)
    nLCUTS::Int                         # number of Lazy Cuts (calls of function)
    nUCUTS::Int                         # number of User Cuts (calls of function)
    nLCUTSC::Int                        # number of Lazy Cuts (cplex output)
    nUCUTSC::Int                        # number of User Cuts (cplex output)
    nLITR::Int                          # number of iterations of Lazy CB (without abort)
    nUITR::Int                          # number of iterations of User CB (without abort)
    nUCBaborts::Int                     # number of aborts in user callback
    nBB::Int                            # number of Branch & Bound nodes
    nSingletons::Int                    # number of nodes which are a singleton in at least one scenario
    rtCPLEX::Float64                    # runtime: cplex
    rtRootCuts::Float64                 # runtime: manually solve root node before B&B
    rtBM::Float64                       # runtime: build cplex model
    rtLG::Float64                       # runtime: create live-graphs
    rtTMIN::Float64                     # runtime: calculating TLmin
    rtR::Float64                        # runtime: calculating set R
    rtMNC::Float64                      # runtime: calculating set marginal node contributions
    rtCSOL::Float64                     # runtime: check solution
    rtLazyCB::Float64                   # runtime: lazy cuts
    rtUserCB::Float64                   # runtime: user cuts
    rtTOT::Float64                      # runtime: total
    rtEVAL::Float64                     # runtime: evaluating the solutions on larger set of scenarios
    rtEVALMAR::Float64                  # runtime: evaluating the solutions on larger set of scenarios
    rtHMAR::Float64                     # runtime: marginal gain heuristic
    rtHBIN::Float64                     # runtime: best individual heuristic
    rtHbc::Float64                      # runtime: betwenness centrality heuristic
    rtHcc::Float64                      # runtime: closeness centrality heuristic
    rtHpr::Float64                      # runtime: reverse page rank heuristic
    rtHtr::Float64                      # runtime: best tunk rank heuristic
    rtHdc::Float64                      # runtime: best degree centrality heuristic
    rtHed::Float64                      # runtime: best expected degree heuristic
    rtHrm::Float64                      # runtime: best retweets heuristic
    rtHmx::Float64                      # runtime: matrix powering heuristic
    memMaxUse::Int                      # mem: maximum usage
    memLG::Int                          # mem: consumption of live-graphs (tuple)
    memTLmin::Int                       # mem: TLmin matrix
    memR::Int                           # mem: inst.valnh BitArray
    nR::Int                             # number of elements in reachability set
    nLG::Int                            # number of elements in live graphs
    aggoutput::Int                      # aggregated output (0=single iteration of SAA, 1=aggregated output of an SAA)
    mySSAItr::Int                       # number of current SAA iteration
    sol_valid::Int                      # solution valid
    exitflag::Int                       # exit status [1 = safe exit, 2=timeoverflow 3=memoverflow, 4=unbounded, 5=infeasible, 6=no integer solution yet, 7=cpx:out-of-memory -1=unknown state]
    exitflagLP::Int                     # exit status [1 = safe exit, 2=timeoverflow -1=unknown state]
    best_ss_id::Int                     # SSA iteration in which the best seed set was found
    best_ss_idMAR::Int                  # SSA iteration in which the best seed set was found
    nexceptions::Int                    # number of exceptions during run
    abs_rank_bc::Vector{Tuple{Int,Int}} # absolute ranks of seed nodes according to betweenness_centrality
    abs_rank_dc::Vector{Tuple{Int,Int}} # absolute ranks of seed nodes according to outdegree_centrality
    abs_rank_cc::Vector{Tuple{Int,Int}} # absolute ranks of seed nodes according to closeness_centrality
    abs_rank_ed::Vector{Tuple{Int,Int}} # absolute ranks of seed nodes according to expected outdegree
    abs_rank_pr::Vector{Tuple{Int,Int}} # absolute ranks of seed nodes according to pagerank
    abs_rank_rm::Vector{Tuple{Int,Int}} # absolute ranks of seed nodes according to retweet/mentions (as in http://www.webecologyproject.org/wp-content/uploads/2009/09/influence-report-final.pdf)
    abs_rank_tr::Vector{Tuple{Int,Int}} # absolute ranks of seed nodes according to tunkRank
    abs_rank_bcMAR::Vector{Tuple{Int,Int}}      # absolute ranks of seed nodes according to betweenness_centrality
    abs_rank_dcMAR::Vector{Tuple{Int,Int}}      # absolute ranks of seed nodes according to outdegree_centrality
    abs_rank_ccMAR::Vector{Tuple{Int,Int}}      # absolute ranks of seed nodes according to closeness_centrality
    abs_rank_edMAR::Vector{Tuple{Int,Int}}      # absolute ranks of seed nodes according to expected outdegree
    abs_rank_prMAR::Vector{Tuple{Int,Int}}      # absolute ranks of seed nodes according to page rank
    abs_rank_rmMAR::Vector{Tuple{Int,Int}}      # absolute ranks of seed nodes according to replies and mentions
    abs_rank_trMAR::Vector{Tuple{Int,Int}}      # absolute ranks of seed nodes according to tunk rank
    rel_rank_bc::Vector{Tuple{Int,Float64}}     # relative ranks of seed nodes according to betweenness_centrality
    rel_rank_dc::Vector{Tuple{Int,Float64}}     # relative ranks of seed nodes according to outdegree centrality
    rel_rank_cc::Vector{Tuple{Int,Float64}}     # relative ranks of seed nodes according to closeness_centrality
    rel_rank_ed::Vector{Tuple{Int,Float64}}     # relative ranks of seed nodes according to expected outdegree
    rel_rank_pr::Vector{Tuple{Int,Float64}}     # relative ranks of seed nodes according to reverse page rank
    rel_rank_rm::Vector{Tuple{Int,Float64}}     # relative ranks of seed nodes according to retweets and mentions
    rel_rank_tr::Vector{Tuple{Int,Float64}}     # relative ranks of seed nodes according to tunk rank
    rel_rank_bcMAR::Vector{Tuple{Int,Float64}}  # relative ranks of seed nodes according to betweenness_centrality
    rel_rank_dcMAR::Vector{Tuple{Int,Float64}}  # relative ranks of seed nodes according to outdegree centrality
    rel_rank_ccMAR::Vector{Tuple{Int,Float64}}  # relative ranks of seed nodes according to closeness_centrality
    rel_rank_edMAR::Vector{Tuple{Int,Float64}}  # relative ranks of seed nodes according to expected outdegree
    rel_rank_prMAR::Vector{Tuple{Int,Float64}}  # relative ranks of seed nodes according to reverse page rank
    rel_rank_rmMAR::Vector{Tuple{Int,Float64}}  # relative ranks of seed nodes according to retweets and mentions
    rel_rank_trMAR::Vector{Tuple{Int,Float64}}  # relative ranks of seed nodes according to tunk rank
    stage::String                       # stage
    msg::String                         # exit message (e.g., error description)
end


# CREATE EMPTY RESULTS
function createEmptyResults( params::Dict, inst::instance )
    # create 'empty' results
    res::results = results( params["opath"],
                   params["ifile"],
                   "-", "-",
                   zeros(Float64, params["nSSAItr"]),
                   zeros(Float64, params["nSSAItr"]),
                   -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0,
                   -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0,
                   -ones(Float64, 3), -ones(Float64, 3),
                   0, 0, 0, 0, 0, 0, 0, 0, -1, 0,
                   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                   0, 0, 0, 0, 0, 0,
                   0, 1, -1, -1, -1, -1, -1, 0,
                   Tuple{Int,Int}[], Tuple{Int,Int}[], Tuple{Int,Int}[], Tuple{Int,Int}[], Tuple{Int,Int}[], Tuple{Int,Int}[], Tuple{Int,Int}[],
                   Tuple{Int,Int}[], Tuple{Int,Int}[], Tuple{Int,Int}[], Tuple{Int,Int}[], Tuple{Int,Int}[], Tuple{Int,Int}[], Tuple{Int,Int}[],
                   Tuple{Int,Float64}[], Tuple{Int,Float64}[], Tuple{Int,Float64}[], Tuple{Int,Float64}[], Tuple{Int,Float64}[], Tuple{Int,Float64}[], Tuple{Int,Float64}[],
                   Tuple{Int,Float64}[], Tuple{Int,Float64}[], Tuple{Int,Float64}[], Tuple{Int,Float64}[], Tuple{Int,Float64}[], Tuple{Int,Float64}[], Tuple{Int,Float64}[],
                   "-","-")
    # set output filename
    res.ofile = ( params["ifile"] *
                  "-nsc-"  * string( inst.nsc ) *
                  "-kF-"   * string( inst.kF  ) *
                  "-kL-"   * string( inst.kL  ) *
#                  "-ssl-"  * string( params["ssltype"]  ) *
                  "-m-"    * string( inst.mtype  ) *
#                  "-benf-" * string( params["benf"]  ) *
                  "-pp-"   * string( params["prepro"]  ) *
#                  "-hcb-"  * string( params["hcallback"]  ) *
                  "-ht-"   * string( params["htype"]  ) *
                  "-sr-"   * string( params["solveRoot"]  ) *
                  "-rc-"   * string( params["rootCuts"]  ) *
                  "-rca-"   * string( params["rootCutsAlpha"]  ) *
                  "-rcl-"   * string( params["rootCutsLambda"]  ) *
                  "-ca-"   * string( params["cutType"]  ) *
                  "-lu-"   * string( params["leaderUtility"]  ) *
#                  "-oa-"   * string( params["OAtype"]  ) *
#                  "-lpm-"   * string( params["LPMethod"]  ) *
#                  "-ca-"   * string( params["cutAlphaType"]  ) *
#                  "-cb-"   * string( params["cutBetaType"]  ) *
                  "-fr-"   * string( params["fractionOfResistantNodes"]  ) *
                  "-pV-"   * string( params["probV"]  ) *
                  "-ct-"   * string( params["relcuttol"]  ) *
                  "-ctU-"   * string( params["relcuttolU"]  ) )
    return res::results
end


# BFS DATA STRUCT
mutable struct bfsdata
    cur_level_L::Vector{Tnodes}         #current level leader
    cur_level_F::Vector{Tnodes}         #current level follower
    next_level_L::Vector{Tnodes}        #current level leader
    next_level_F::Vector{Tnodes}        #current level follower
    tvisited::Vector{UInt8}             #activated by follower(=1)/leader(=2)
    bvisited::Vector{Bool}              #node visited or not
    xv::Vector{Float64}                 # count viewings (possibly fractional)
    xv_val::Matrix{Float64}             # values xv_val[i,w] retrieved from callback for OA
    xf::Vector{Float64}                 # count forwards (possibly fractional)
    xfb::Vector{Bool}                   # forwarders or not (binary used for integer solutions)
    xin::Vector{Float64}                # interior point for in-out method
    yv::Vector{Float64}                 # count viewings (possibly fractional)
    yv_val::Matrix{Float64}             # values yv_val[i,w] retrieved from callback for OA
    yf::Vector{Float64}                 # count forwards (possibly fractional)
    yfb::Vector{Bool}                   # forwarders or not (binary used for integer solutions)
    alpha::Vector{Float64}              # Lagrange multiplier
    phi::Vector{Float64}                # Lagrange multiplier
    pi::Vector{Float64}                 # Lagrange multiplier
    rho::Vector{Float64}                # Benders cut coefficients
    Phi::Float64                        # current solution of benders subproblem
end


# CREATE EMPTY BFS DATA STRUCT
function createEmptyBFSdata( nnodes::Int, nsc::Int )
    bfs::bfsdata = bfsdata( zeros(Tnodes, nnodes),
                    zeros(Tnodes, nnodes),
                    zeros(Tnodes, nnodes),
                    zeros(Tnodes, nnodes),
                    zeros(UInt8,  nnodes),
                    zeros(Bool,   nnodes),
                    zeros(Float64, nnodes),
                    zeros(Float64, nnodes, nsc),
                    zeros(Float64, nnodes),
                    zeros(Bool, nnodes),
                    zeros(Float64, nnodes),
                    zeros(Float64, nnodes),
                    zeros(Float64, nnodes, nsc),
                    zeros(Float64, nnodes),
                    zeros(Bool, nnodes),
                    zeros(Float64, nnodes),
                    zeros(Float64, nnodes),
                    zeros(Float64, nnodes),
                    zeros(Float64, nnodes),
                    0.0,
                     )
    return bfs::bfsdata
end

# CLEAR BFS DATA
function clear!( b::bfsdata )
    empty!(b.cur_level_L)
    empty!(b.cur_level_F)
    empty!(b.next_level_L)
    empty!(b.next_level_F)
    fill!(b.tvisited, zero(UInt8))
    fill!(b.bvisited, zero(Bool))
    fill!(b.xv, zero(Float64))
    fill!(b.xf, zero(Float64))
    fill!(b.xfb, zero(Bool))
    fill!(b.yv, zero(Float64))
    fill!(b.yf, zero(Float64))
    fill!(b.yfb, zero(Bool))
    fill!(b.alpha, zero(Float64))
    fill!(b.phi, zero(Float64))
    fill!(b.pi, zero(Float64))
    fill!(b.rho, zero(Float64))
    b.Phi = 0
    return nothing
end


#RESET RESULTS AFTER COMPUTING LEADERS SEED SET
function reset_results!( res::results )
        res.exitflag    = -1
        res.msg         = "-"
        res.cpxStatus   = "-"
        res.sol_valid   = -1
        res.LPrelax     = -1
        res.bestbound   = -1
        res.obj_F       = -1
        res.obj_L       = -1
        res.gap         = -1
        res.nBENCUTS    = 0
        res.nLCUTS      = 0
        res.nUCUTS      = 0
        res.nBB         = 0
        res.nSingletons = 0
        res.rtTMIN      = 0
        res.rtR         = 0
        res.rtMNC       = 0
        res.rtBM        = 0
        res.rtCPLEX     = 0
        res.rtLazyCB    = 0
        res.rtUserCB    = 0
        res.memMaxUse   = 0
        res.memR        = 0
        return nothing
end


#RESET RESULTS AFTER COMPUTING LEADERS SEED SET
function reset_results_exceptions!( res::results )
        res.exitflag    = -1
        res.sol_valid   = -1
        res.LPrelax     = -1
        res.bestbound   = -1
        res.obj_F       = -1
        res.obj_L       = -1
        res.gap         = -1
        res.nBENCUTS    = 0
        res.nLCUTS      = 0
        res.nUCUTS      = 0
        res.nBB         = 0
        res.nSingletons = 0
        res.rtTMIN      = 0
        res.rtR         = 0
        res.rtMNC       = 0
        res.rtBM        = 0
        res.rtCPLEX     = 0
        res.rtLazyCB    = 0
        res.rtUserCB    = 0
        res.memMaxUse   = 0
        res.memR        = 0
        return nothing
end


