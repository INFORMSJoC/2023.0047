"""
USED PACKAGES: DataFrames CSV Query Colors LaTeXStrings PyPlot Seaborn Pycall

NEEDED INPUT: DataFrame with dummy-all-ones-column named :ones
DataFrame must contain the runtime in a column named :rt

NEEDS PARAMETER STRUCT THAT CONTAINS:
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

NEEDS DICTIONARY THAT MAPS AXIS LABELS, E.G.,
mydict = Dict( :ifile => "Instance",
               :nsc   => "|Omega|",
               :LPgap => "LPgap")

AUTHOR: M.Kahr
"""

#INCLUDES
include("declarations.jl")


#***************************************************************
#********* PERFORMANCE PLOTS ***********************************
#***************************************************************

#PYPLOT: PLOT RUNTIME PERFORMANCE PLOTS GROUPED BY PARAMETER s (e.g., kF, ifile,...)
#DATAFRAME MYDF NEEDS A DUMMY COLUMN :ones (i.e., an all-ones vector)
#option filter is for the case that the data is filterd before and juste extends the filename
function plotRuntimeBySymbol( mydf::DataFrame, p::parameters, s::Symbol, filter::String="")
    ofile = "pp_" * string( s ) * "_rt" * filter
    sort!(mydf, (s))
    s_vals = unique(mydf[!,s]) #param values (e.g., instance names)
    ninst = size(mydf,1)
    groupsize = ninst / length( s_vals )
    mydf.rt = mydf.rt ./ 1000           

    #empty plot
    xlim_lb = 0
    p.pptype == "log" ? xlim_lb = 1 : nothing                       #change lower bound if logarithmic scale
    fig, ax = plt.subplots(nrows=1, ncols=1 ,figsize=(5,5))
#    plt.xlim(xlim_lb, p.tlim+100)
    plt.xlim(xlim_lb, 7.3)
    plt.ylim(0,105)
    ax.set_xscale( p.pptype ) #argument is a string


    #create runtime plots
    for i=1:length(s_vals)
        #create subdf
        df = mydf[(mydf[!,:opt] .== 1) .&
                  (mydf[!,s] .== s_vals[i] ), :]
        sort!(df, (:rt))

        #append dummy column (for ploting lines to the timelimit)
        tempdf = similar(df,1)
#        tempdf[:,:rt]   .= p.tlim
        tempdf[:,:rt]   .= 7.2
        tempdf[:,:ones] .= 0
        append!( df, tempdf )

        #add line to rtplot
        x=df[!,:rt]
        y=(cumsum(df[!,:ones]) ./groupsize) .*100
        plt.plot(x, y,      dashes   = p.dstyle[i],
                            lw       = p.lwidth,
                            label    = string(s_vals[i]),
                            drawstyle="steps-post",
#                            label  = string(s_vals[i]),
#                            smooth = false
                            )
    end

    legendpos = p.legendpos
    #legendpos = "above"

#    plt.xlabel("runtime [s]")
    plt.xlabel("runtime [1000s]")
    plt.ylabel("\\# instances solved [%]")
    yname = ""
    haskey( p.mydict, s ) ? yname = p.mydict[s] : xcat_name = string( s )
    if legendpos == "above"
        plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc="lower left", ncol=3, fontsize="small", mode="expand", borderaxespad=0.)
    else
        plt.legend(title=yname, loc="lower right",fancybox="false", framealpha=0.9)
    end
#    plt.legend(title=yname, fancybox="false", framealpha=0.9, bbox_to_anchor=(1.04,1))
#    plt.tight_layout(rect=[0,0,0.75,1])
#    plt.tight_layout() # avoid cutting of labels
    plt.tick_params(right=true)
#    plt.savefig(ofile * ".pdf")
    plt.savefig(p.opath * ofile)
    display(fig)
    p.saveTikz ? tikzplotlib.save(p.opath * ofile * ".tex") : nothing
    ax=nothing
    plt.close(fig)



    # gap plots
    ofile = "pp_" * string( s ) * "_gap" * filter
    fig = plt.figure(ofile ,figsize=(5,5))
#    plt.xlim(0, 1)
    plt.ylim(0,105)
    plt.xlim(0, 3)
    plt.xticks(np.arange(0, 3, step=0.5))
    if s == :mtype 
        plt.xlim(0, 105)
        plt.xticks(np.arange(0, 105, step=10))
    end
    #plt.xlim(0, 105)
    #plt.xticks(np.arange(0, 105, step=10))
    for i=1:length(s_vals)
        #create subdf
#        df = mydf[ (mydf[!,:mtype] .== 5) .& (mydf[!,s] .== s_vals[i] ), :]
        df = mydf[ (mydf[!,s] .== s_vals[i] ), :]
        sort!(df, (:gap))

        #cumsum to the first nonzero gap value (=offset). Remark: needed that plots don't start from zero
        csum = 0
        for l=1:size( df, 1 )
            if( df[!,:opt][l] == 1 )
                df[!,:ones][l] = 0
                csum += 1
            end
        end

        #add line to gapplot
        x=df[!,:gap]
        y=(( cumsum(df[!,:ones]) .+ csum) ./groupsize) .*100
        plt.plot(x, y,      dashes   = p.dstyle[i],
                            lw       = p.lwidth,
                            label    = string(s_vals[i]),
                            drawstyle="steps-post",
                            )
    end

    #add plot attributes
    plt.xlabel("optimality gap [%]")
    plt.ylabel("\\#instances solved [%]")
    if legendpos == "above"
        plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc="lower left", ncol=3, fontsize="small", mode="expand", borderaxespad=0.)
    else
        plt.legend(title=yname, loc="lower right",fancybox="false", framealpha=0.9)
    end
    plt.savefig(p.opath * ofile)
    p.saveTikz ? tikzplotlib.save(p.opath * ofile * ".tex") : nothing
    display(fig)
    plt.close(fig)

    println("called print function")
    return nothing
end


#data=value of interest, s=group of interest
function cumulativePlotsBySymbols( mydf::DataFrame, p::parameters, data::Symbol, s::Symbol, filter::String="")
    ofile = "pp_" * string( data ) * "_by_" * string( s ) * filter
    sort!(mydf, (s))
    s_vals = unique(mydf[!,s]) #(e.g., instance names)
    maxlength = 0
    for subdf in groupby(mydf, s)
        if nrow(subdf) > maxlength
            maxlength = nrow(subdf)
        end
    end
    println(maxlength)

    #empty plot
    xlim_lb = 0
    p.pptype == "log" ? xlim_lb = 1 : nothing                       #change lower bound if logarithmic scale
    fig, ax = plt.subplots(nrows=1, ncols=1 ,figsize=(5,5))
    plt.xlim(xlim_lb, 105)
    plt.ylim(0,105)
    ax.set_xscale( p.pptype ) #argument is a string


    #create runtime plots
    for i=1:length(s_vals)
        #create subdf
        df = mydf[mydf[!,s] .== s_vals[i], :]
        sort!(df, (:obj_ratio))

        #add line to rtplot
        x=(collect(1:nrow(df)) ./ maxlength) .* 100
        y=df[!,data]

        plt.plot(x, y,      dashes   = p.dstyle[i],
                            lw       = p.lwidth,
                            label    = string(s_vals[i]),
                            drawstyle="steps-post",
                            )
    end

    #legendpos = "lower right"
    legendpos = "above"
    if p.legendpos == "above"
#        plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc="lower left", ncol=3, fontsize="small", mode="expand", borderaxespad=0.)
    elseif p.legendpos == "right"
        plt.legend(bbox_to_anchor=(1.02, 0.9, 0.8, .102), loc=2, mode="expand", borderaxespad=0.,frameon=false)
    else
#        plt.legend(title=yname, loc="lower right",fancybox="false", framealpha=0.9)
    end


    plt.xlabel("\\# instances solved [%]")
#    plt.xlabel("runtime [1000s]")
    plt.ylabel("Heuristic solution quality [%]")
    yname = ""
    haskey( p.mydict, s ) ? yname = p.mydict[s] : xcat_name = string( s )
    plt.tick_params(right=true)
#    plt.savefig(ofile * ".pdf")
    plt.savefig(p.opath * ofile)
    display(fig)
    p.saveTikz ? tikzplotlib.save(p.opath * ofile * ".tex") : nothing
    ax=nothing
    plt.close(fig)

    println("called print function")
    return nothing
end



#*********************************************************
#********* BOX PLOTS *************************************
#*********************************************************

#PYPLOT: GENERIC BOXPLOTS
function simpleBoxplot( rdf::DataFrame, p::parameters, xcat::Symbol, ycat::Symbol )
    ofile = "bp_" * string( ycat ) * "_by_" * string( xcat )
    sort!(rdf, (xcat) )

    ordc2 = ["F","O","T", "R25", "R50", "R75"]                                   #NOTE: problem specific
    orderdict = Dict(x => i for (i,x) in enumerate(ordc2))
    sort!(rdf; cols = [xcat], by = x->orderdict[x])

    #PROBLEM specific settings: ...........rename categories if this is more meaningful
    xcat == :gdensity ? xcat = :ifile : nothing
    xcat == :narcs    ? xcat = :ifile : nothing
    xcat == :aDeg     ? xcat = :ifile : nothing
    xcat == :aEDeg    ? xcat = :ifile : nothing
    xcat == :aRsize   ? xcat = :ifile : nothing
    #..........................................

    # get categories
    s_vals = unique(rdf[!,xcat]) #param ausprägungen (e.g., instance names)
    data = Array{Any}(undef, length(s_vals))
    mylabels = String[]

    for i=1:length(s_vals)
        data[i] = rdf[ rdf[!,xcat] .== s_vals[i], :][!,ycat]
        push!( mylabels, string(s_vals[i]) )
    end

    fig = plt.figure(ofile ,figsize=(5,5))
    plt.boxplot(data, labels = mylabels, showfliers=p.showfliers)
#    plt.boxplot(data, labels = mylabels)
    xcat_name = ""
    ycat_name = ""
    haskey( p.mydict, xcat ) ? xcat_name = p.mydict[xcat] : xcat_name = string( xcat )
    haskey( p.mydict, ycat ) ? ycat_name = p.mydict[ycat] : ycat_name = string( ycat )
    #p.crun == 2       ? plt.ylim(0,1) : nothing                                 #TODO: Problem specific!!
    plt.xlabel( xcat_name )
    plt.ylabel( ycat_name )
    xcat == :ifile ? plt.xticks( rotation = 90 ) : nothing
    plt.tight_layout()                                              #avoid cutting of instance names
    plt.savefig(p.opath * ofile)
    p.saveTikz ? tikzplotlib.save(p.opath * ofile * ".tex") : nothing
    display(fig)
    plt.close(fig)
    return nothing
end



#*********************************************************
#********* SWARM PLOTS ***********************************
#*********************************************************

#PYPLOT / SEABORN: GENERIC SWARMPLOT
function simpleSwarmPlot( mydf::AbstractDataFrame, p::parameters, xcat::Symbol, ycat::Symbol, ext::String="", myhue::Symbol=:dummy )
    ofile = "swp_" * string( ycat ) * "_by_" * string( xcat ) * ext
    sort!(mydf, (xcat) )
    xcat_vals = unique(mydf[!,xcat])   #var ausprägungen (e.g., instance names)

    fig = plt.figure(ofile ,figsize=(5,5))
    if myhue == :dummy
        sns.swarmplot(x=mydf[!,xcat], y=mydf[!,ycat], color="#8192a1")
    else
        sns.swarmplot(x=mydf[!,xcat], y=mydf[!,ycat], hue = mydf[!,myhue], dodge = p.mydodge )
    end
    xcat_name = ""
    ycat_name = ""
    haskey( p.mydict, xcat ) ? xcat_name = p.mydict[xcat] : xcat_name = string( xcat )
    haskey( p.mydict, ycat ) ? ycat_name = p.mydict[ycat] : ycat_name = string( ycat )
    plt.xlabel( xcat_name )
    plt.ylabel( ycat_name )
    plt.ylim(15, 105)
    if p.legendpos == "above"
        plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc="lower left", ncol=3, fontsize="small", mode="expand", borderaxespad=0.)
    end
    if p.legendpos == "right"
        plt.legend(bbox_to_anchor=(1.02, 0.9, 0.8, .102), loc=2, mode="expand", borderaxespad=0.,frameon=false)
#        plt.legend().set_visible(false)
    end
    xcat == :ifile ? plt.xticks( rotation = 90 ) : nothing
    plt.tight_layout()                                              #avoid cutting of instance names
    plt.savefig(p.opath * ofile)
    p.saveTikz ? tikzplotlib.save(p.opath * ofile * ".tex") : nothing
    display(fig)
    plt.close(fig)
end


function simpleSwarmBoxPlot( mydf::AbstractDataFrame, p::parameters, xcat::Symbol, ycat::Symbol, m::String="", myhue::Symbol=:dummy )
    ofile = "swbp_" * string( ycat ) * "_by_" * string( xcat ) * "_by_" * m
    sort!(mydf, (xcat) )

    #TODO: this produces error in ANAL_IMP
    if m ∉ ["R", "LU"]
        ordc2 = ["F","O","T", "R25", "R50", "R75"]                                   #NOTE: problem specific
        orderdict = Dict(x => i for (i,x) in enumerate(ordc2))
        sort!(mydf; cols = [xcat], by = x->orderdict[x])
    end

    xcat_vals = unique(mydf[!,xcat])   #var ausprägungen (e.g., instance names)

    fig = plt.figure(ofile ,figsize=(5,5))
    if m == "R"
        plt.ylim(-0.2, 3.5)
    elseif m == "LU"
        plt.ylim(-5, 105)
    else
        plt.ylim(-5, 300)
    end
    if myhue == :dummy
        sns.boxplot(x=mydf[!,xcat], y=mydf[!,ycat], palette="Blues")
        sns.swarmplot(x=mydf[!,xcat], y=mydf[!,ycat], color="Black")
    else
        sns.swarmplot(x=mydf[!,xcat], y=mydf[!,ycat], hue = mydf[!,myhue], dodge = p.mydodge )
    end
    xcat_name = ""
    ycat_name = ""
    haskey( p.mydict, xcat ) ? xcat_name = p.mydict[xcat] : xcat_name = string( xcat )
    haskey( p.mydict, ycat ) ? ycat_name = p.mydict[ycat] : ycat_name = string( ycat )
    plt.xlabel( xcat_name )
    plt.ylabel( ycat_name )
    if p.legendpos == "above"
        plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc="lower left", ncol=3, fontsize="small", mode="expand", borderaxespad=0.)
    end
    xcat == :ifile ? plt.xticks( rotation = 90 ) : nothing
    plt.tight_layout()                                              #avoid cutting of instance names
    plt.savefig(p.opath * ofile * ".pdf")
    p.saveTikz ? tikzplotlib.save(p.opath * ofile * ".tex") : nothing
    display(fig)
    plt.close(fig)
end
