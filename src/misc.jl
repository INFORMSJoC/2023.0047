# MISC FUNCTIONS
#include("var_declarations.jl")


function setInstancePathAndName!( params::Dict )
#    params["ifile"] ="soc-advogato.im"
#    params["ifile"] ="soc-anybeat.im"
#    params["ifile"] ="msg-college.im"
#    params["ifile"] ="msg-email-eu.im"
#    params["ifile"] ="TW-austria-20200814-Y2020-anonymized"
#    params["ifile"] ="TW-giftideas-20200810-Y2020-anonymized"
#    params["ifile"] ="TW-greenenergy-20200811-Y2020-anonymized"
#    params["ifile"] ="TW-organicfood-20200813-Y2020-anonymized"
#    params["ifile"] ="TW-orms-20200807-Y2020-anonymized"
#    params["ifile"] ="TW-skateboarding-20200810-Y2020-anonymized"
    params["ifile"] = "TW-travelling-20200807-Y2020-anonymized"
#    params["ifile"] = "D-n25-k4-b0.3-beta2.0-18.0-i1"
end

# GET MEMORY CONSUMPTION
function get_mem_use()
    f::IOStream         = open( "/proc/self/stat", "r" )
    s::AbstractString   = read( f, String )
    vsize::Int          = parse( Int64, split( s )[23] )
    mb::Int             = Int( ceil( vsize / ( 1024 * 1024 ) ) )
    close(f)
    return mb::Int
end

# MEMORY CONSUMPTION TEST
function memOK( mb::Int, memlimit::Int )
    if( mb > memlimit )
        println("### Memory-usage too high: $mb (MB)")
        return false
    else
        return true
    end
end


# COMPILER HELP (e.g., cast x::Any to x::Float64)
function annotate( a::Array{Float64,1} )
    return a::Array{Float64,1}
end

function annotateMatrix( a::Array{Float64,2} )
    return a::Array{Float64,2}
end

function stringify( s::Any )
    return s::String
end

function floatify( s::Any )
    return s::Float64
end


# HANDLE EXCEPTIONS
function handleExceptions( inst::instance, res::results, params::Dict, wo::Bool, wherefrom::String, msg::Any )
    res.nexceptions += 1
    reset_results_exceptions!( res )
    if( typeof(msg) == CPLEX.CplexError )       #cplex error
        if msg.code == 1001                                                 #cplex hit physical limit
            res.exitflag = 7
            res.msg = "Code 1001: cpx out-of-memory"
            wo ? writeOutputToCSV( params, inst, res ) : nothing
            println("### CPLEX ERROR")
#            exit()
        end
        if msg.code == 1811                                                 #cplex hit physical limit
            res.exitflag = 7
            res.msg = "Code 1811: Attempt to invoke unsupported operation"
            wo ? writeOutputToCSV( params, inst, res ) : nothing
            println("### CPLEX ERROR")
#            exit()
        end
    elseif( typeof(msg) == OutOfMemoryError )   #OutOfMemoryError()
            res.exitflag = 8
            res.msg = string(wherefrom, ": OutOfMemoryError", )
            wo ? writeOutputToCSV( params, inst, res ) : nothing
            println("### OUT OF MEMORY ERROR")
#            exit()
    else                                        #some other error
            strmsg = string( msg )
            res.exitflag = -1
            res.msg = string( wherefrom, ": ", strmsg )
            wo ? writeOutputToCSV( params, inst, res ) : nothing
            println("### SOME UNKNOWN ERROR")
#            exit()
    end
    println("### UNKNOWN ERROR IN:", wherefrom )
    return nothing
end

