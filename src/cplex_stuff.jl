#****************************************************************
#**** CPLEX RELATED STUFF
#****************************************************************
#DEFINE CPLEX CONSTANTS
const CPXCALLBACKINFO_THREADID  = convert(Cint,0)
const CPXCALLBACKINFO_NODECOUNT = convert(Cint,1)
const CPXCALLBACKINFO_ITCOUNT   = convert(Cint,2)
const CPXCALLBACKINFO_BEST_SOL  = convert(Cint,3)
const CPXCALLBACKINFO_BEST_BND  = convert(Cint,4)
const CPXCALLBACKINFO_THREADS   = convert(Cint,5)
const CPXCALLBACKINFO_FEASIBLE  = convert(Cint,6)
const CPXCALLBACKINFO_TIME      = convert(Cint,7)
const CPXCALLBACKINFO_DETTIME   = convert(Cint,8)


#GETTERS
#current b&b node
function cpx_get_nodecount( cb_data::CPLEX.CallbackContext )
    data_p = Ref{Cint}()
    return_status = CPLEX.@cpx_ccall(callbackgetinfoint,
                                        Cint,
                                        (Ptr{Cvoid}, Cint, Ref{Cint}),
                                        cb_data.context, CPXCALLBACKINFO_NODECOUNT, data_p)
    cur_node = data_p[]

    if return_status != 0
        @warn "error retrieving current B&B node"
    end
    return cur_node::Int32
end

#current current best bound
function cpx_get_best_bnd( cb_data::CPLEX.CallbackContext )
    data_p = Ref{Cdouble}()
    return_status = CPLEX.@cpx_ccall(callbackgetinfodbl,
                                        Cint,
                                        (Ptr{Cvoid}, Cint, Ref{Cdouble}),
                                        cb_data.context, CPXCALLBACKINFO_BEST_BND, data_p)
    best_bnd = data_p[]

    if return_status != 0
        @warn "error retrieving best_bound"
    end
    return best_bnd::Float64
end

#current current best sol
function cpx_get_best_sol( cb_data::CPLEX.CallbackContext )
    data_p = Ref{Cdouble}()
    return_status = CPLEX.@cpx_ccall(callbackgetinfodbl,
                                        Cint,
                                        (Ptr{Cvoid}, Cint, Ref{Cdouble}),
                                        cb_data.context, CPXCALLBACKINFO_BEST_SOL, data_p)
    best_sol = data_p[]

    if return_status != 0
        @warn "error retrieving best_sol"
    end
    return best_sol::Float64
end

#compute current mip gap [%]
function cpx_get_mip_gap( cb_data::CPLEX.CallbackContext )
    mip_gap::Float64 = 0.0
    best_bnd = cpx_get_best_bnd( cb_data )
    best_sol = cpx_get_best_sol( cb_data )
    mip_gap  = ((best_bnd - best_sol) / (0.0000000001 + best_sol)) * 100
    return mip_gap::Float64
end

#get number of cuts
function cpx_get_num_cuts( model::Model, cuttype::Int32 )
    num_cuts::Int32 = -1
    moi_model       = backend( model )                     # MOI
    cplex_model     = moi_model.inner                      # CPLEX.jl
    num_cuts        = CPLEX.get_num_cuts( cplex_model, cuttype )
    return num_cuts::Int32
end


# CREATE AN EMPTY MODEL (JuMP 0.21)
function createEmptyModelMOI( params::Dict, timelimit::Int )
    #create model
    model = JuMP.direct_model( CPLEX.Optimizer() )

    #set params
    if params["solveRoot"]  #deactivate all heuristics etc. to obtain the true LP bound (Not implemented here)
        #general
        set_optimizer_attribute( model, "CPXPARAM_TimeLimit", timelimit )
        set_optimizer_attribute( model, "CPXPARAM_WorkMem", params["memlimit"] )
        set_optimizer_attribute( model, "CPXPARAM_Threads", 1 )
        set_optimizer_attribute( model, "CPXPARAM_MIP_Display", 2 )
        #set_optimizer_attribute(model, "CPXPARAM_MIP_Limits_Nodes", 0 )

        #numerics
        set_optimizer_attribute( model, "CPXPARAM_MIP_Tolerances_MIPGap", 0.00001 )
        set_optimizer_attribute( model, "CPXPARAM_Emphasis_Numerical", 1 )

        #used LP methods (0:auto (default), 1:primal simplex, 2:dual simplex, 4:barrier)
        set_optimizer_attribute( model, "CPXPARAM_LPMethod", params["LPMethod"] )
        set_optimizer_attribute( model, "CPXPARAM_MIP_Strategy_SubAlgorithm", params["LPMethod"] )
        set_optimizer_attribute( model, "CPXPARAM_MIP_Strategy_StartAlgorithm", params["LPMethod"] )

        #deactivate general purpose cuts
        set_optimizer_attribute( model, "CPXPARAM_MIP_Limits_EachCutLimit", 0 )
        set_optimizer_attribute( model, "CPXPARAM_MIP_Cuts_Gomory", -1 )
        set_optimizer_attribute( model, "CPXPARAM_MIP_Cuts_LiftProj", -1 )
        set_optimizer_attribute( model, "CPXPARAM_MIP_Strategy_CallbackReducedLP", 0 )

        #deactivate preprocessing
        set_optimizer_attribute( model, "CPXPARAM_Preprocessing_Presolve", 0 )
        set_optimizer_attribute( model, "CPXPARAM_Preprocessing_Relax", 0 )
        set_optimizer_attribute( model, "CPXPARAM_Preprocessing_RepeatPresolve", 0 )
        set_optimizer_attribute( model, "CPXPARAM_MIP_Strategy_PresolveNode", -1 )
        set_optimizer_attribute( model, "CPXPARAM_MIP_Strategy_Probe", -1 )

        #deactivate  heuristics
        set_optimizer_attribute( model, "CPXPARAM_MIP_Strategy_HeuristicFreq", -1 )
        set_optimizer_attribute( model, "CPXPARAM_MIP_Strategy_RINSHeur", -1 )
        set_optimizer_attribute( model, "CPXPARAM_MIP_Strategy_FPHeur", -1 )
        set_optimizer_attribute( model, "CPXPARAM_MIP_Strategy_LBHeur", 0 )
    else
        #general
        set_optimizer_attribute( model, "CPXPARAM_TimeLimit", timelimit )
        set_optimizer_attribute( model, "CPXPARAM_WorkMem", params["memlimit"] )
        set_optimizer_attribute( model, "CPXPARAM_Threads", 1 )
        set_optimizer_attribute( model, "CPXPARAM_MIP_Display", 2 )
        #set_optimizer_attribute(model, "CPXPARAM_MIP_Limits_Nodes", 0 )

        #numerics
        set_optimizer_attribute( model, "CPXPARAM_MIP_Tolerances_MIPGap", 0.00001 )
        #set_optimizer_attribute( model, "CPXPARAM_Emphasis_Numerical", 1 )

        #used LP methods (0:auto (default), 1:primal simplex, 2:dual simplex, 4:barrier)
        #set_optimizer_attribute( model, "CPXPARAM_LPMethod", params["LPMethod"] )
        #set_optimizer_attribute( model, "CPXPARAM_MIP_Strategy_SubAlgorithm", params["LPMethod"] )
        #set_optimizer_attribute( model, "CPXPARAM_MIP_Strategy_StartAlgorithm", params["LPMethod"] )
    end

    return model::Model
end


#****************************************************************
#**** CALLBACK FOR CHECKING MEMORY CONSUMPTION DURING SOLVING 
#****************************************************************
# PROCSTATUS - check if the used memory is with in the set limit, utherwise abort solving (packed in heuristic callback)
function startProcStatus( cb_data::CPLEX.CallbackContext, inst::instance, memorylimit::Int, res::results )
    mb = get_mem_use()
    if( mb > res.memMaxUse )
        res.memMaxUse = mb
        #println("MaxMem = ", mb)
    end
    if( !memOK(mb, memorylimit) )
        res.exitflag = 3
        return CPLEX.terminate( cb_data.model )
    end
    return nothing
end