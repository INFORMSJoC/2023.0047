# PARSE PARAMETER INPUTS
using ArgParse
#include("var_declarations.jl")

# returns a dictionary with problem parameters
function parse_commandline()
    s = ArgParseSettings(description = "Program description",
                         version="0.1")

    @add_arg_table! s begin
        "--root_path"
            help = "root folder: /my/path/to/root_folder/"
            arg_type = String
            default = "/home/kahr/svn/projects/misc/logit_code_public/"
        "--ipath", "-p"
            help = "path to instance file: /my/path/to/instance/"
            arg_type = String
            #default = "/home/kahr/svn/projects/SocialNetworks/instances/smallworld-ic/"
            default = "/home/kahr/svn/projects/misc/logit_code_public/data/instances/"
        "--opath"
            help = "path to result folder: /my/path/to/results/"
            arg_type = String
            default = "/home/kahr/svn/projects/misc/logit_code_public/results/"
        "--ifile", "-f"
            help = "instance file name, e.g.,: TW-orms-20200807-Y2020-anonymized"
            arg_type = String
            default = "D-n25-k4-b0.3-beta2.0-18.0-i1"
        "--itype"
            help = "instance type: 1=artificial, 2=SNAP, 3=Twitter"
            arg_type = Int64
            default = 3
        "--mtype"
            help = "model type (1 = max. forwarders (variant F), 2=max. total views (variant T), 3=max. organic reach (variant O), 4=max. patronage with outer approximation (variant R), 5=max. patronage with generalized Benders decomposition (variant R)"
            arg_type = Int64
            default = 5
        "--benf"
            help = "if true: add benders fractional cuts"
            arg_type = Bool
            default = true
        "--nSSAItr"
            help = "number of sample average iterations"
            arg_type = Int64
            default = 2
        "--Nscenarios"
            help = "number of scenarios (evaluation) |Ω''|"
            arg_type = Int64
            default = 100
        "--nscenarios", "-w"
            help = "number of scenarios |Ω'|"
            arg_type = Int64
            default = 10
        "--kL"
            help = "seed-set size of L(leader). Note: only seed sets for variant R with |L| ∈ {5,10,15} are precomputed. Otherwise a seed set using the Local-Index-Rank heuristic is computed (cf., Liu2017)"
            arg_type = Int
            default = 0
        "--kF"
            help = "seed-set size of F(ollower)"
            arg_type = Int64
            default = 5
        "--cutType"
            help = "unused: Benders cut type from [1:9] (see sheet of paper) (5 is std setting)"
            arg_type = Int
            default = 5
        "--prepro"
            help = "preprocessing true or false"
            arg_type = Bool
            default = true
        "--solveRoot"
            help = "unused: deactivate CPLEX heuristics etc. to get the true LP bound"
            arg_type = Bool
            default = false
        "--rootManualLPrelax"
            help = "unused: manually solve LP relaxation"
            arg_type = Bool
            default = false
        "--rootCuts"
            help = "unused: add Benders cuts in root node manually"
            arg_type = Bool
            default = false
        "--rootCutsAlpha"
            help = "unused: α ∈ [0,1] defining the convex combination of vectors for in-out method"
            arg_type = Float64
            default = 0.3
        "--rootCutsLambda"
            help = "unused: λ ∈ [0,1] defining the convex combination of vectors for in-out method"
            arg_type = Float64
            default = 0.5
        "--debug", "-d"
            help = "debug (or verbose) level from 0 to 5."
            arg_type = Int
            default = 1
        "--console", "-c"
            help = "set to false if program is started externally from script (or shell); set to true if started from REPL"
            arg_type = Bool
            default = true
        "--precompile"
            help = "if true: start a dummy instance to get everything compiled, then start the real instance"
            arg_type = Bool
            default = false
        "--writeoutput"
            help = "if true: write result files"
            arg_type = Bool
            default = true
        "--resistanceHurdle"
            help = "fractions below this value will never patronize"
            arg_type = Float64
            default = 0.10
        "--fractionOfResistantNodes"
            help = "fraction fraction of resistant nodes"
            arg_type = Float64
            default = 0.25
        "--relcuttol"
            help = "relative tolerance for adding lazy Benders cuts (+1)"
            arg_type = Float64
            default = 1.005
        "--relcuttolU"
            help = "relative tolerance for adding user Benders cuts (+1)"
            arg_type = Float64
            default = 1.01
        "--probF"
            help = "probability for activation arcs (if not sampled)"
            arg_type = Float64
            default = 0.1
        "--leaderUtility"
            help = "factor of scaling leader utility"
            arg_type = Int
            default = 1
        "--probV"
            help = "probability for passivizing arcs (if not sampled)"
            arg_type = Float64
            default = 0.05
        "--ssltype"
            help = "unused: seed-set type of leader"
            arg_type = Int
            default = 4
        "--computeSL"
            help = "unused: write leader seedset to file"
            arg_type = Bool
            default = false
        "--timelimit"
            help = "timelimit for CPLEX (seconds)"
            arg_type = Int
            default = 10
        "--memlimit"
            help = "memory limit (MB)"
            arg_type = Int
            default = 12000
        "--rseed"
            help = "random seed"
            arg_type = Int
            default = 1
        "--ptype"
            help = "unused"
            arg_type = Int
            default = 1
        "--htype"
            help = "heuristic type: 0=no, 3=MAR"
            arg_type = Int
            default = 3
        "--hcallback"
            help = "unused: if true: use heuristic callback"
            arg_type = Bool
            default = false
        "--OAtype"
            help = "unused: 0=∀ω ∈ Ω in RHS, 1=sum_{ω ∈ Ω} in RHS (only one var in obj function)"
            arg_type = Int
            default = 0
        "--LPMethod"
            help = "unused: Method CPLEX uses for solving LPs; 0:auto (default), 1:primal simplex, 2:dual simplex, 4:barrier"
            arg_type = Int
            default = 0
        "--solveLP"
            help = "unused: solve only LP relaxation (true/false)"
            arg_type = Bool
            default = false
        "--exportGraphML"
            help = "unused: export to graphML (true/false)"
            arg_type = Bool
            default = false
        "--computeInstanceStats"
            help = "unused: compute instance stats (e.g., avg deg,..)"
            arg_type = Bool
            default = false
        "--computeRsize"
            help = "unused: compute instance stats (e.g., avg deg,..)"
            arg_type = Bool
            default = false
        "--clevel"
            help = "unused: confidence level of solution quality"
            arg_type = Float64
            default = 0.1
        "--computeRankings"
            help = "compute rankings of seed sets (e.g., in betweenness centrality)"
            arg_type = Bool
            default = true
    end

    return parse_args(s)
end

#function main()
#    params = parse_commandline()
#    println("Parsed args:")
#    for (arg,val) in params
#        println("  $arg  =>  $val")
#    end
#end

#main()
