
# LOAD PACKAGES
using DataFrames
using CSV
using Query
using Statistics
using StatsBase
using Distributions
using Colors
using LaTeXStrings
import PyPlot
import Seaborn
using PyCall

#set paths
path  = ""
file  = ""
#set outpath
opath = joinpath(dirname(@__FILE__), "plots/")
isdir(opath) ? nothing : mkdir(opath)

#import matplotlib.pyplot
plt = pyimport("matplotlib.pyplot")
sns = pyimport("seaborn")
np  = pyimport("numpy")
tikzplotlib = pyimport("tikzplotlib")

#set matplotlib params
rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
rcParams["font.family"]    = "serif"
rcParams["savefig.format"] = "png"
rcParams["hist.bins"]      = 100
rcParams["legend.framealpha"] = 0.9


#DEFINE PARAMETER STRUCT
mutable struct parameters
    tlim::Int64                 #time limit
    opath::String               #outpath
    ninst::Int64                #number of instances
    inames::Array{String,1}     #vector of instance names
    colors::Array{Color,1}      #vector of colors
    dstyle::Vector{Vector{Int}} #vector of line styles (used in "dashes" keyword argument)
    lwidth::Any                 #line thickness
    mydict::Dict{Symbol,String} #converts a column name to a string
    fontfamily::String          #font type of figure labels
    granPP::Int64               #granularity performance plot
    granGAP::Float64            #granularity gap plots
    saveTikz::Bool              #save plot as tikz or not
    showfliers::Bool            #show outliers in boxplot or not
    pptype::String              #string that indicates perfanceplot type "linear", "log", "logit"
    legendpos::String           #position of the legend ("best", "above", "lower left", "upper left"...) -> should be extended
    mydodge::Bool               #whether or not data gets dodged in swarmplot (i.e., similar data is placed vertically)
    xlimLB::Int                 # x-axis limit lower bound
    xlimUB::Int                 # x-axis limit upper bound
    ylimLB::Int                 # y-axis limit lower bound
    ylimUB::Int                 # y-axis limit upper bound
end


#DEFINE LEGENDS
mydict = Dict( :ifile           => "Instance",
               :nsc             => L"$|\Omega'|$",
               :kF              => L"$|F|$",
               :kL              => L"$|L|$",
               :obj_F           => L"$\sigma(L,F)$",
               :obj_L_eval      => L"$\lambda_{\Omega''}^{\mathrm{BEN}}(L,F)$",
               :obj_L           => L"$\lambda_{\Omega''}^{\mathrm{BEN}}(L,F)$",
               :obj_L_evalChen  => L"$\lambda_{\Omega''}^{\mathrm{MAR}}(L,F)$",
               :obj_L_evalIND   => L"$\lambda_{\Omega''}^{\mathrm{BIN}}(L,F)$",
               :obj             => L"$\sigma_{\Omega''}(L,F)$",
               :obj_F_evalChen  => L"$\sigma_{\Omega''}^{\mathrm{MAR}}(L,F)$",
               :obj_F_evalIND   => L"$\sigma_{\Omega''}^{\mathrm{BIN}}(L,F)$",
               :obj_total_F         => L"$\sigma_{\Omega''}^{\mathrm{T}}(L,F)$",
               :obj_organic_F       => L"$\sigma_{\Omega''}^{\mathrm{O}}(L,F)$",
               :obj_resistant_F25   => L"$\sigma_{\Omega''}^{\mathrm{R25}}(L,F)$",
               :obj_resistant_F50   => L"$\sigma_{\Omega''}^{\mathrm{R50}}(L,F)$",
               :obj_resistant_F75   => L"$\sigma_{\Omega''}^{\mathrm{R75}}(L,F)$",
               :loss_total_F         => L"$\sigma_{\Omega''}^{\mathrm{T}}(L,F)$",
               :loss_organic_F       => L"$\sigma_{\Omega''}^{\mathrm{O}}(L,F)$",
               :loss_forwards_F       => L"$\sigma_{\Omega''}^{\mathrm{F}}(L,F)$",
               :loss_resistant_F25   => L"$\sigma_{\Omega''}^{\mathrm{R25}}(L,F)$",
               :loss_resistant_F50   => L"$\sigma_{\Omega''}^{\mathrm{R50}}(L,F)$",
               :loss_resistant_F75   => L"$\sigma_{\Omega''}^{\mathrm{R75}}(L,F)$",
               :setSimOverModels    => L"$\cap_{i} \frac{F_i}{|F|} [\%]$",
               :gap_approx_rel  => L"$\Delta [\%]$",
               :rt              => L"$\bar{rt}^{\mathrm{BEN}}$",
               :rtHChen         => L"$\bar{rt}^{\mathrm{MAR}}$",
               :rtHIND          => L"$\bar{rt}^{\mathrm{BIN}}$",
               :artHChen        => L"$\bar{rt}^{\mathrm{MAR}}$",
               :artHIND         => L"$\bar{rt}^{\mathrm{BIN}}$",
               :mtypeString     => "Metric",
               :cdl90 => "confidence level 90%",
               :cdl95 => "confidence level 95%",
               :cdl99 => "confidence level 99%",
               :setSim => L"$\cap_{i} F_i / |F| [\%]$",
               :best_ss_id => L"\# SSA iterations for $F^*$",
               :EDeg  => L"Expected node degrees $\mathbb{|\delta(i)|}$",
               :Deg  =>  L"Node degrees $|\delta(i)|$",
               :relL  => L"$1- \hat{\sigma}_{\Omega''}(L, F^*) / \hat{\sigma}_{\Omega''}(L, \emptyset)$",
               :relF  => L"$1- \hat{\sigma}_{\Omega''}(L, F^*) / \hat{\sigma}_{\Omega''}(\emptyset, F^*)$",
#               :rt    => "runtime [s]",
               :htype => "heuristic type",
    #           :obj_F_eval => "GB",
               :obj_Fbc_eval => "BC",
               :obj_Fdc_eval => "DC",
               :obj_Fed_eval => "ED",
               :obj_FMAR_eval => "MG",
               :obj_Fpr_eval => "PR",
               :obj_Frm_eval => "RM",
               :obj_Ftr_eval => "TR",
    #           :rtCPLEX => "GB",
               :rtHbc => "BC",
               :rtHdc => "DC",
               :rtHed => "ED",
               :rtHMAR => "MG",
               :rtHpr => "PR",
               :rtHrm => "RM",
               :rtHtr => "TR",
               :obj_F_eval => L"$e^{-\infty}$",
               :obj_F_eval_1 => L"$e^0$",
               :obj_F_eval_2 => L"$e^1$",
               :obj_F_eval_3 => L"$e^2$",
               :obj_F_eval_4 => L"$e^3$",
               :obj_F_eval_5 => L"$e^4$",
               :rtCPLEX => L"$e^{-\infty}$",
               :rtCPLEX_1 => L"$e^0$",
               :rtCPLEX_2 => L"$e^1$",
               :rtCPLEX_3 => L"$e^2$",
               :rtCPLEX_4 => L"$e^3$",
               :rtCPLEX_5 => L"$e^4$",
               :nBB   => "#BB nodes",
               :LPgap => "LPgap")


#DEFINE COLORS
colors = ["blue", "darkred", "darkgreen", "snow4", "darkorange3", "brown", "purple2", "green", "gold4", "grey", "olive", "violetred1"]
colvec = Color[]
for i=1:length(colors) #convert colors for color blind people
    col = parse( Colorant, colors[i] )

    colblindcol = protanopic( col )
    push!( colvec, colblindcol )
end

#define empty variables
xparams          = Symbol[]
xgroups          = Symbol[]
yparams          = Symbol[]
paramsOfInterest = Symbol[]

#RENAME INSTANCES
function rename_instances!( rdf::DataFrame )
    for l=1:size(rdf,1)
        rdf.ifile[l] == "TW-datascience-20190414-Y2019"   ? rdf.ifile[l] = "tw-datascience"   : nothing
        rdf.ifile[l] == "TW-giftideas-20191010-Y2019"     ? rdf.ifile[l] = "tw-giftideas"     : nothing
        rdf.ifile[l] == "TW-nrw2019-20191009-Y2019"       ? rdf.ifile[l] = "tw-nrw2019"       : nothing
        rdf.ifile[l] == "TW-orms-20191009-Y2019"          ? rdf.ifile[l] = "tw-orms"          : nothing
        rdf.ifile[l] == "TW-valentinesday-20191010-Y2019" ? rdf.ifile[l] = "tw-valentinesday" : nothing
        rdf.ifile[l] == "TW-vienna-20190412-Y2019"        ? rdf.ifile[l] = "tw-vienna"        : nothing
        rdf.ifile[l] == "TW-ootd-20191012-Y2019"          ? rdf.ifile[l] = "tw-ootd"          : nothing
        rdf.ifile[l] == "msg-college.im"                  ? rdf.ifile[l] = "msg-college"      : nothing
        rdf.ifile[l] == "msg-email-eu.im"                 ? rdf.ifile[l] = "msg-email-eu"     : nothing
        rdf.ifile[l] == "p2p-gnutella09.im"               ? rdf.ifile[l] = "p2p-gnutella09"   : nothing
        rdf.ifile[l] == "p2p-gnutella30.im"               ? rdf.ifile[l] = "p2p-gnutella30"   : nothing
        rdf.ifile[l] == "soc-advogato.im"                 ? rdf.ifile[l] = "soc-advogato"     : nothing
        rdf.ifile[l] == "soc-anybeat.im"                  ? rdf.ifile[l] = "soc-anybeat"      : nothing
        rdf.ifile[l] == "TW-travelling-20200807-Y2020"    ? rdf.ifile[l] = "tw-travelling"    : nothing
        rdf.ifile[l] == "TW-orms-20200807-Y2020"          ? rdf.ifile[l] = "tw-orms"          : nothing
        rdf.ifile[l] == "TW-skateboarding-20200810-Y2020" ? rdf.ifile[l] = "tw-skateboarding" : nothing
        rdf.ifile[l] == "TW-naturelovers-20200809-Y2020"  ? rdf.ifile[l] = "tw-naturelovers"  : nothing
        rdf.ifile[l] == "TW-giftideas-20200810-Y2020"     ? rdf.ifile[l] = "tw-giftideas"     : nothing
        rdf.ifile[l] == "TW-greenenergy-20200811-Y2020"   ? rdf.ifile[l] = "tw-greenenergy"   : nothing
        rdf.ifile[l] == "TW-organicfood-20200813-Y2020"   ? rdf.ifile[l] = "tw-organicfood"   : nothing
        rdf.ifile[l] == "TW-austria-20200814-Y2020"       ? rdf.ifile[l] = "tw-austria"       : nothing
    end
end


function addInstanceStats!( rdf::DataFrame )
    #nodes
    statspath = "/home/kahr/svn/projects/SocialNetworks/codeComp/results/inst_stats/"
    statsfile = "all-inststatsnode_twprob.csv"

    #read and aggregate data
    statsdf = CSV.read( statspath * statsfile; copycols=true, delim=';' )
    rename_instances!( statsdf )
    aggdf = by(statsdf, :ifile, :Deg => mean, :EDeg => mean)

    #add columns
    rdf.aDeg = zeros( Float64, nrow(rdf) )                          #average node degree
    rdf.aEDeg = zeros( Float64, nrow(rdf) )                         #average expected node degree

    iname::String = ""
    for l=1:size(rdf,1)
        iname       = rdf.ifile[l]
        rdf.aDeg[l] = aggdf[(aggdf.ifile .== iname), :Deg_mean][1]
        rdf.aEDeg[l]= aggdf[(aggdf.ifile .== iname), :EDeg_mean][1]
    end


    #edges
    statspath = "/home/kahr/svn/projects/SocialNetworks/codeComp/results/inst_stats/"
    statsfile = "all-inststatsedge_twprob.csv"

    #read and aggregate data
    statsdf = CSV.read( statspath * statsfile; copycols=true, delim=';' )
    rename_instances!( statsdf )
    aggdf = by(statsdf, :ifile, :p => mean, :p => sum)

    #add columns
    rdf.ap   = zeros( Float64, nrow(rdf) )                          #average edge probability
    rdf.sump = zeros( Float64, nrow(rdf) )                          #average edge probability

    iname = ""
    for l=1:size(rdf,1)
        iname       = rdf.ifile[l]
        rdf.ap[l]   = aggdf[(aggdf.ifile .== iname), :p_mean][1]
        rdf.sump[l] = aggdf[(aggdf.ifile .== iname), :p_sum][1]
    end
end
