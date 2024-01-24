# WRITE RESULTS TO CSV FILE
function writeOutputToCSV( params::Dict, inst::instance, res::results )
      # write headlines file
      hname = "0_headlines.csv"
      f = open( joinpath(res.opath, hname), "w")
          print(f,
                "ifile;",                             # instance filename
                "exitflag;",                          # exitflag ( 1 = everything ok )
                "exitflagLP;",                        # exitflag LP relaxation ( 1 = everything ok )
                "stage;",                             # last stage of the program (for debugging)
                "msg;",                               # optional message
                "aggoutput;",                         # aggregated output file (0=false: output of one SAA iteration, 1=true: aggregated file with average values of SAA iterations)
                "mySSAItr;",                          # current SAA iteration (for debugging)
                "tlim;",                              # set time limit
                "nsc;",                               # |Ω|
                "Nsc;",                               # |Ω'| evaluation scenarios
                "nSSAItr;",                           # number of set SAA iterations
                "kF;",                                # |S^F| cardinality of follower seed set
                "kL;",                                # |S^L| cardinality of leader seed set
                "nnodes;",                            # |V|
                "narcs;",                             # |A|
#                "avg_deg;",                           # average node degree
                "hashval;",                           # hash value describing instance
                "probV;",                             # viewing probability
                "cutType;",                           # Benders cut types
#                "cutA;",                              # Benders cut α types
#                "cutB;",                              # Benders cut β types
                "mtype;",                             # model type
                "HCB;",                               # use heuristic callback
                "htype;",                             # type of heuristic
                "BENF;",                              # BENFractional cuts
                "PREPRO;",                            # preprocessing true/false
                "OAtype;",                            # OA type (forall vs sum)
                "LPMethod;",                          # LPMethod used by CPLEX
                "solveRoot;",                         # deactivate CPLEX heuristics etc. to get true LP bound
                "rootCuts;",                          # manually solve root node as LP
                "rootCutsAlpha;",                     # parameter α for manual separation in root node
                "rootCutsLambda;",                    # parameter λ for manual separation in root node
                "resistanceHurdle;",                  # selfexplaining
                "fractionOfResistantNodes;",          # selfexplaining
                "cpxStatus;",                         # cplex status
                "cpxStatusLP;",                       # cplex status LP relaxation
                "sol_valid;",                         # 1=solution valid, 0=otherwise
                "LPrelax;",                           # LP relaxation
                "LPgap;",                             # LP relaxation gap [%]
                "bestbound;",                         # best objective value bound
                "obj_F;",                             # objective value (of F(ollower))
                "obj_L;",                             # objective value of Leader - obtained from propagation
                "obj_F_solcheck;",                    # objective obtained from solution checker
                "obj_F_eval;",                        # best objective value of follwer after evalution
                "obj_L_eval;",                        # best objective value of follwer after evalution
                "fcut_error;",                        # relative error between solution checker and cplex (in percent)
                "gap;",                               # optimality gap (in percent)
                "gap_approx_abs;",                    # estimated absolute approximation gap
                "gap_approx_rel;",                    # estimated relative approximation gap (in percent)
                "best_ss_id;",                        # id of the SSA iteration in which best seed set was found
                "nBENCUTS;",                          # number of benders cuts (function calls)
                "nLCUTS;",                            # number of lazy benders cuts (function calls)
                "nUCUTS;",                            # number of user benders cuts (function calls)
                "nLCUTSC;",                           # number of lazy benders cuts (cplex output excluding purged cuts)
                "nUCUTSC;",                           # number of user benders cuts (cplex output excluding purged cuts)
                "nUCBaborts;",                         # number of aborts in UCB
                "nBB;",                               # number of B&Bnodes
                "nSingletons;",                       # number of nodes which are a singeltons in at least one scenario
                "obj_FMAR;",                          # objective value MAR heuristic
                "obj_FMAR_eval;",                     # objective value MAR heuristic
                "obj_LMAR_eval;",                     # objective value MAR heuristic
                "obj_Fbc_eval;",                      # objective value MAR heuristic
                "obj_Fcc_eval;",                      # objective value MAR heuristic
                "obj_Fdc_eval;",                      # objective value MAR heuristic
                "obj_Fed_eval;",                      # objective value MAR heuristic
                "obj_Fpr_eval;",                      # objective value MAR heuristic
                "obj_Frm_eval;",                      # objective value MAR heuristic
                "obj_Ftr_eval;",                      # objective value MAR heuristic
                "obj_Lbc_eval;",                      # objective value MAR heuristic
                "obj_Lcc_eval;",                      # objective value MAR heuristic
                "obj_Ldc_eval;",                      # objective value MAR heuristic
                "obj_Led_eval;",                      # objective value MAR heuristic
                "obj_Lpr_eval;",                      # objective value MAR heuristic
                "obj_Lrm_eval;",                      # objective value MAR heuristic
                "obj_Ltr_eval;",                      # objective value tunkrank
                "rtLG;",                              # rt: create live graphs ## ALL rt measured in seconds
#                "rtTMAX;",                            # rt: get Tmax(ω)_i, ∀ i,ω,companies
                "rtTMIN;",                            # rt: get TLmin
                "rtR;",                               # rt: calculate valid Neiborhood
                "rtMNC;",                             # rt: marginal node contribution
                "rtBM;",                              # rt: build model
                "rtCPLEX;",                           # rt: solve cplex
                "rtRootCuts;",                        # rt: manually solving root node before B&B
                "rtLazyCB;",                          # rt: lazy cuts cplex
                "rtUserCB;",                          # rt: user cuts cplex
                "rtEVAL;",                            # rt: function evaluations
                "rtHMAR;",                            # rt: heuristic eval
                "rtHbc;",                               # rt: heuristic eval
                "rtHcc;",                               # rt: heuristic eval
                "rtHdc;",                               # rt: heuristic eval
                "rtHed;",                               # rt: heuristic eval
                "rtHpr;",                               # rt: heuristic eval
                "rtHrm;",                               # rt: heuristic eval
                "rtHtr;",                               # rt: heuristic eval
#                "rtCSOL[s];",                         # rt: check solution
#                "rtTOT[s];",                          # rt: total
                "memMaxUse;",                         # mem: maximum usage ## ALL Measured in MiB
                "memLG;",                             # mem: livegraphs struct (tuples)
                "memR;",                              # mem: valid R^ω BitArray
                "nLG;",                              # number of elements in live-graphs
                "nR;",                               # number of elements reachability sets
                "nexceptions;",                       # mem: valid R^ω BitArray
                "relcuttol;",                         # relative fractional cut tolerance (+1)
                "relcuttolU;",                         # relative fractional cut tolerance (+1)
                "abs_rank_bc;",                         #absolute rank betweenness
                "abs_rank_cc;",                         #absolute rank closeness
                "abs_rank_dc;",                         #absolute rank outdegree
                "abs_rank_ed;",                         #absolute rank expected outdegree
                "abs_rank_pr;",                         #absolute rank expected outdegree
                "abs_rank_rm;",                         #absolute rank expected outdegree
                "abs_rank_tr;",                         #absolute rank expected outdegree
                "abs_rank_bcMAR;",                      #absolute rank betweenness
                "abs_rank_ccMAR;",                      #absolute rank closeness
                "abs_rank_dcMAR;",                      #absolute rank outdegree
                "abs_rank_edMAR;",                      #absolute rank expected outdegree
                "abs_rank_prMAR;",                      #absolute rank expected outdegree
                "abs_rank_rmMAR;",                      #absolute rank expected outdegree
                "abs_rank_trMAR;",                      #absolute rank expected outdegree
                "rel_rank_bc;",                         #relative rank betweenness
                "rel_rank_cc;",                         #relative rank closeness
                "rel_rank_dc;",                         #relative rank outdegree
                "rel_rank_ed;",                         #relative rank expected outdegree
                "rel_rank_pr;",                         #relative rank expected outdegree
                "rel_rank_rm;",                         #relative rank expected outdegree
                "rel_rank_tr;",                         #relative rank expected outdegree
                "rel_rank_bcMAR;",                      #relative rank betweenness
                "rel_rank_ccMAR;",                      #relative rank closeness
                "rel_rank_dcMAR;",                      #relative rank outdegree
                "rel_rank_edMAR;",                      #relative rank expected outdegree
                "rel_rank_prMAR;",                      #relative rank expected outdegree
                "rel_rank_rmMAR;",                      #relative rank expected outdegree
                "rel_rank_trMAR;",                      #relative rank expected outdegree
                "SL;",                                # seed set leader
                "SFMAR;",                             # seed set MAR heuristic
                "SFbc;",                              # seed set MAR heuristic
                "SFcc;",                             # seed set MAR heuristic
                "SFdc;",                             # seed set MAR heuristic
                "SFed;",                             # seed set MAR heuristic
                "SFpr;",                             # seed set MAR heuristic
                "SFrm;",                             # seed set MAR heuristic
                "SFtr;",                             # seed set MAR heuristic
                "tnameL;",                            # twitter names Leader
                "tnameF;",                            # twitter names Follower
#                "nTotalViewsF;",                      # twitter names Follower
#                "nOrganicViewsF;",                    # twitter names Follower
#                "nForwardsF;",                       # twitter names Follower
#                "nTotalViewsL;",                     # twitter names Follower
#                "nOrganicViewsL;",                   # twitter names Follower
#                "nForwardsL;",                       # twitter names Follower
                "SF;",                                # seed set follower
                #"SF\n"                                # seed set follower
                "leaderUtility\n"                      # seed set follower
                )
      close(f)

    # write results to csv
    mySSAItr = "-"
    if( res.mySSAItr < inst.nSSAItr )
          mySSAItr = "0" * string(res.mySSAItr)
    else
          mySSAItr = string(res.mySSAItr)
    end

    tnamesL = String[]
    tnamesF = String[]
    if(params["itype"] == 3)
          tnamesL = inst.tname[inst.SL]
          tnamesF = inst.tname[inst.SF]
    end

    out_filename = res.ofile * "-SSAItr-" * mySSAItr * "-agg-" * string(res.aggoutput)
    f = open( joinpath(res.opath, out_filename), "w")
        print(f,
              inst.ifile, ";",                      # instance filename
              res.exitflag, ";",                    # exitflag ( 1 = everything ok )
              res.exitflagLP, ";",                  # exitflag ( 1 = everything ok )
              res.stage, ";",                       # last stage of the program
              res.msg, ";",                         # error message
              res.aggoutput, ";",                   # aggregated output (0=false, 1=true)
              res.mySSAItr, ";",                    # current iteration of SSA
              params["timelimit"], ";",             # time limit
              inst.nsc, ";",                        # |Ω|
              inst.Nsc, ";",                        # |Ω'| evaluation scenarios
              inst.nSSAItr, ";",                    # number of SSA Iterations
              inst.kF, ";",                         # |S^F|
              params["kL"], ";",                    # |S^L|
              inst.nnodes, ";",                     # |V|
              inst.narcs, ";",                      # |A|
#              inst.avg_deg, ";",                    # average node degree
              inst.hashval, ";",                    # hashvalue describing instance
              params["probV"], ";",                 # viewing probability
              params["cutType"], ";",               # Benders cut α type
#              params["cutAlphaType"], ";",          # Benders cut α type
#              params["cutBetaType"], ";",           # Benders cut β type
              inst.mtype, ";",                      # model type
              params["hcallback"], ";",             # use heuristic callback (true/false)
              params["htype"], ";",                 # heuristic type
              params["benf"], ";",                  # include benders fractional cuts (true/false)
              params["prepro"], ";",                # preprocessing (true/false)
              params["OAtype"], ";",                # OAtype (forall vs. sum RHS)
              params["LPMethod"], ";",              # LPMethod used by CPLEX
              params["solveRoot"], ";",             # deactivate CPLEX heursitics etc. to get true LP bound
              params["rootCuts"], ";",              # manually solve root node as LP
              params["rootCutsAlpha"], ";",         # parameter α ∈ [0,1] for manual separation in root node
              params["rootCutsLambda"], ";",        # parameter λ ∈ [0,1] for manual separation in root node
              params["resistanceHurdle"], ";",      # self explaining
              params["fractionOfResistantNodes"], ";", # self explaining
              res.cpxStatus, ";",                   # cplex status
              res.cpxStatusLP, ";",                 # cplex status LP relaxation
              res.sol_valid, ";",                   # solution valid
              res.LPrelax, ";",                     # LP relaxation
              res.LPgap, ";",                       # LP gap
              res.bestbound, ";",                   # best objective value bound
              res.obj_F, ";",                       # obj_F val (model)
              res.obj_L, ";",                       # objective value of Leader - obtained from propagation
              res.obj_F_solcheck, ";",              # objective value obtained from solution check
              res.obj_F_eval, ";",                  # best objective value of follower after evaluation
              res.obj_L_eval, ";",                  # objective value of leader at evaluation
              res.fcut_error, ";",                  # relative error between solution of checker and cplex
              res.gap, ";",                         # optimality gap (cplex)
              res.gap_approx_abs, ";",              # estimated absolute approximation gap
              res.gap_approx_rel, ";",              # estimated relative approximation gap
              res.best_ss_id, ";",                  # id of best seed set in SSA iteration
              res.nBENCUTS, ";",                    # number of benders cuts (function calls)
              res.nLCUTS, ";",                      # number of lazy benders cuts (function calls)
              res.nUCUTS, ";",                      # number of user benders cuts (function calls)
              res.nLCUTSC, ";",                     # number of lazy benders cuts (cplex output excluding purged cuts)
              res.nUCUTSC, ";",                     # number of user benders cuts (cplex output excluding purged cuts)
              res.nUCBaborts, ";",                  # number of aborts in ucb
              res.nBB, ";",                         # number of B&B nodes
              res.nSingletons, ";",                 # number of nodes that are a singleton in at least one scenario
              res.obj_FMAR, ";",                    # objective value MAR heuristic
              res.obj_FMAR_eval, ";",               # objective value MAR heuristic
              res.obj_LMAR_eval, ";",               # objective value MAR heuristic
              res.obj_Fbc_eval, ";",                # objective value bteweeness heuristic
              res.obj_Fcc_eval, ";",                # objective value closeness heuristic
              res.obj_Fdc_eval, ";",                # objective value degree heuristic
              res.obj_Fed_eval, ";",                # objective value expected outdegree heuristic
              res.obj_Fpr_eval, ";",                # objective value expected outdegree heuristic
              res.obj_Frm_eval, ";",                # objective value expected outdegree heuristic
              res.obj_Ftr_eval, ";",                # objective value expected outdegree heuristic
              res.obj_Lbc_eval, ";",                # objective value bteweeness heuristic
              res.obj_Lcc_eval, ";",                # objective value closeness heuristic
              res.obj_Ldc_eval, ";",                # objective value degree heuristic
              res.obj_Led_eval, ";",                # objective value expected outdegree heuristic
              res.obj_Lpr_eval, ";",                # objective value expected outdegree heuristic
              res.obj_Lrm_eval, ";",                # objective value expected outdegree heuristic
              res.obj_Ltr_eval, ";",                # objective value expected outdegree heuristic
              res.rtLG, ";",                        # rt: create live graphs
#              res.rtTMAX, ";",                      # rt: get Tmax(ω)_i, ∀ i,ω,companies
              res.rtTMIN, ";",                      # rt: get TLmin
              res.rtR, ";",                         # rt: get valid Neighborhood R^ω
              res.rtMNC, ";",                       # rt: marginal node contributions
              res.rtBM, ";",                        # rt: build model
              res.rtCPLEX, ";",                     # rt: solve cplex
              res.rtRootCuts, ";",                  # rt: manually solve root node before B&B
              res.rtLazyCB, ";",                    # rt: lazy CB
              res.rtUserCB, ";",                    # rt: User CB
              res.rtEVAL, ";",                      # rt: Evaluation
              res.rtHMAR, ";",                      # rt: heristic MAR
              res.rtHbc, ";",                         # rt: heristic MAR
              res.rtHcc, ";",                         # rt: heristic MAR
              res.rtHdc, ";",                         # rt: heristic MAR
              res.rtHed, ";",                         # rt: heristic MAR
              res.rtHpr, ";",                         # rt: heristic MAR
              res.rtHrm, ";",                         # rt: heristic MAR
              res.rtHtr, ";",                         # rt: heristic MAR
#              res.rtCSOL, ";",                      # rt: check solution
#              res.rtTOT, ";",                       # rt: total
              res.memMaxUse, ";",                   # mem: maximum usage
              res.memLG, ";",                       # mem: live graphs struct (tuples)
              res.memR, ";",                        # mem: valid N^ω BitArray
              res.nLG, ";",                         # number of elements in LG
              res.nR, ";",                          # number of elements in R
              res.nexceptions, ";",                 # mem: number of exceptions in SSA
              params["relcuttol"], ";",             # relative cut tolerance
              params["relcuttolU"], ";",             # relative fractional cut tolerance
              res.abs_rank_bc, ";",                 # absolute rank betwennness
              res.abs_rank_cc, ";",                 # absolute rank closeness
              res.abs_rank_dc, ";",                 # absolute rank degree
              res.abs_rank_ed, ";",                 # absolute rank expected outdegree
              res.abs_rank_pr, ";",                 # absolute rank expected outdegree
              res.abs_rank_rm, ";",                 # absolute rank expected outdegree
              res.abs_rank_tr, ";",                 # absolute rank expected outdegree
              res.abs_rank_bcMAR, ";",                 # absolute rank betwennness
              res.abs_rank_ccMAR, ";",                 # absolute rank closeness
              res.abs_rank_dcMAR, ";",                 # absolute rank degree
              res.abs_rank_edMAR, ";",                 # absolute rank expected outdegree
              res.abs_rank_prMAR, ";",                 # absolute rank expected outdegree
              res.abs_rank_rmMAR, ";",                 # absolute rank expected outdegree
              res.abs_rank_trMAR, ";",                 # absolute rank expected outdegree
              res.rel_rank_bc, ";",                 # relative rank betwennness
              res.rel_rank_cc, ";",                 # relative rank closeness
              res.rel_rank_dc, ";",                 # relative rank degree
              res.rel_rank_ed, ";",                 # relative rank expected outdegree
              res.rel_rank_pr, ";",                 # relative rank expected outdegree
              res.rel_rank_rm, ";",                 # relative rank expected outdegree
              res.rel_rank_tr, ";",                 # relative rank expected outdegree
              res.rel_rank_bcMAR, ";",                 # relative rank betwennness
              res.rel_rank_ccMAR, ";",                 # relative rank closeness
              res.rel_rank_dcMAR, ";",                 # relative rank degree
              res.rel_rank_edMAR, ";",                 # relative rank expected outdegree
              res.rel_rank_prMAR, ";",                 # relative rank expected outdegree
              res.rel_rank_rmMAR, ";",                 # relative rank expected outdegree
              res.rel_rank_trMAR, ";",                 # relative rank expected outdegree
              inst.SL, ";",                         # seed set leader
              inst.SFMAR, ";",                         # seed set MAR
              inst.SFbc, ";",                         # seed set betweenness_centrality
              inst.SFcc, ";",                         # seed set closeness
              inst.SFdc, ";",                         # seed set degree
              inst.SFed, ";",                         # seed set expected degree
              inst.SFpr, ";",                         # seed set pagerank
              inst.SFrm, ";",                         # seed set retweet/mentions
              inst.SFtr, ";",                         # seed set tunkrank
              tnamesL, ";",                         # twitter names leader
              tnamesF, ";",                         # twitter names follower
#              inst.nTotalViewsF, ";",                         # seed set leader
#              inst.nOrganicViewsF, ";",                         # seed set MAR
#              inst.nForwardsF, ";",                         # seed set betweenness_centrality
#              inst.nTotalViewsL, ";",                         # seed set closeness
#              inst.nOrganicViewsL, ";",                         # seed set degree
#              inst.nForwardsL, ";",                         # seed set degree
              inst.SF, ";",                         # seed set follower
#              inst.SF, "\n"                         # seed set follower
              inst.leaderUtility, "\n"                         # seed set follower
               )
    close(f)
end

#write leaders seed set to file
function writeOutputLeadersSeedSet( params::Dict, inst::instance, res::results )
      mykF = ""
      if inst.kF < 10
            mykF = "0" * string(inst.kF)
      else
            mykF = string(inst.kF)
      end

      # (Benders)
      out_filename = "SSL-LOG_BEN_" * inst.ifile * "_k=" * mykF * "_m=" * string(inst.mtype) * "_r=" * string(params["fractionOfResistantNodes"]) * "_w=" * string(inst.nsc)
      f = open( joinpath(res.opath, out_filename) , "w")
          print( f, inst.kF )
          for i ∈ inst.SF
                print(f, " ", i )
          end
          print( f, "\n" )
      close(f)
end
