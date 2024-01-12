function insertSpaceAfterThreeDigits( n::Int )
    nstring = string( n )
    mynewString = ""
    if length(nstring) == 4
        mynewString = nstring[1] * "\\," * nstring[2] * nstring[3] * nstring[4]
    end
    if length(nstring) == 5
        mynewString = nstring[1] * nstring[2] * "\\," * nstring[3] * nstring[4] * nstring[5]
    end
    if length(nstring) == 6
        mynewString = nstring[1] * nstring[2] * nstring[3] * "\\," * nstring[4] * nstring[5] * nstring[6]
    end

    return mynewString
end


function printLatexTableIMP( df::DataFrame, p::parameters )
    colnames = names( df )
    nmrows1::Int  = nrow( df ) / length( unique( df.ifile ) ) #number of rows per instance
    nmrows2::Int  = 7 #length( unique( df.kF ) ) #number of rows per instance

    #print table
    println("\n")
    println(" Printing latex table ")
    print("\n\n\n")

    #commands for rotating entries and top and bottom struts, cf., https://tex.stackexchange.com/questions/65127/extra-vertical-space-after-hline-causes-a-gap-in-the-right-border-of-an-array


    println("\\newcommand{\\STAB}[1]{\\begin{tabular}{@{}c@{}}#1\\end{tabular}}")
    println("\\newcommand\\Tstrut{\\rule{0pt}{2.6ex}}         % = top strut")
    println("\\newcommand\\Bstrut{\\rule[-0.9ex]{0pt}{0pt}}   % = bottom strut")

    #open table
    println( "\\begin{table}[htbp]")
    println( "\t\\caption{mycaption}")
    println( "\t\\begin{adjustbox}{width=0.95\\columnwidth,center}")
    println( "\t\t\\begin{tabular}{crr|rrrrrrrr|rrrrrrrrr}")

    #write headlines
    print("\t\t\t")
    println("& & & \\multicolumn{8}{c|}{runtimes} & \\multicolumn{8}{c}{objective values} \\\\")
    print("\t\t\t")
    for n in colnames
        n == :ifile ? nothing : print(" & ")
        if haskey( p.mydict, n )
            print( p.mydict[n] )
        else
            print( n )
        end
    end
    print("\\Tstrut \\Bstrut \\\\ ")


    #fill data
    #rows
    for i=1:nrow( df )
        (i-1) % nmrows1 == 0 ? print(" \\hline \n") : nothing                   #new outer multi row

        if( (i-1) % nmrows1 != 0 )
            (i-1) % nmrows2 == 0 ? print("\\cline{2-19} \n")  : print(" \n")               #new inner multi row
        else
            #nothing
        end

        print("\t\t\t")

        #columns
        for j=1:ncol( df )
            if j == 1
                (i-1) % nmrows1 == 0 ? print("\\multirow{$nmrows1}{*}{\\STAB{\\rotatebox[origin=c]{90}{", df.ifile[i], "}}}") : print("")     #new outer multi row
            elseif j == 2
                (i-1) % nmrows2 == 0 ? print("& \\multirow{$nmrows2}{*}{", df.kF[i], "}") : print(" & ")        #new inter multi row
            else
                if typeof(df[i,j]) == String
                    print(" & ", df[i,j] )
                else
                    if Int(df[i,j]) >= 1000
                        convertedNumber = insertSpaceAfterThreeDigits( Int(df[i,j]) )
                        print(" & ", convertedNumber )
                    else
                        if Int(df[i,j]) == -969696
                            print(" & ", "-" )
                        else
                            print(" & ", Int(df[i,j]) )
                        end
                    end

                end
            end
        end
#        (i-1) % nmrows1 == 0 ? print("\\Tstrut") : print("")
#        if( (i-1) % nmrows1 != 0 )
#            (i-1) % nmrows2 == 0 ? print("\\Tstrut ")  : print("")
#        end
        (i-1) % nmrows2 == 0 ? print(" \\hspace{-2.5mm} \\Tstrut ")  : print("")
        print(" \\\\")
    end

    print("\n")
    println( "\t\t\\end{tabular}")
    println( "\t\\label{tab:mylabel}")
    println( "\t\\end{adjustbox}")
    println( "\\end{table}")
    println("\n")
end



function printLatexTableCIMP( df::DataFrame, p::parameters )
    colnames = names( df )
    nmrows1::Int  = nrow( df ) / length( unique( df.ifile ) ) #number of rows per instance
    nmrows2::Int  = 7 #length( unique( df.kF ) ) #number of rows per instance

    #print table
    println("\n")
    println(" Printing latex table ")
    print("\n\n\n")

    #commands for rotating entries and top and bottom struts, cf., https://tex.stackexchange.com/questions/65127/extra-vertical-space-after-hline-causes-a-gap-in-the-right-border-of-an-array


    println("\\newcommand{\\STAB}[1]{\\begin{tabular}{@{}c@{}}#1\\end{tabular}}")
    println("\\newcommand\\Tstrut{\\rule{0pt}{2.6ex}}         % = top strut")
    println("\\newcommand\\Bstrut{\\rule[-0.9ex]{0pt}{0pt}}   % = bottom strut")

    #open table
    println( "\\begin{table}[htbp]")
    println( "\t\\caption{mycaption}")
    println( "\t\\begin{adjustbox}{width=0.95\\columnwidth,center}")
    println( "\t\t\\begin{tabular}{lr|rrrrrr|rrrrrr}")

    #write headlines
    print("\t\t\t")
    println("& & \\multicolumn{6}{c|}{runtimes} & \\multicolumn{6}{c}{objective values} \\\\")
    print("\t\t\t")
    for n in colnames
        n == :ifile ? nothing : print(" & ")
        if haskey( p.mydict, n )
            print( p.mydict[n] )
        else
            print( n )
        end
    end
    print("\\Tstrut \\Bstrut \\\\ ")


    #fill data
    #rows
    for i=1:nrow( df )
        (i-1) % nmrows1 == 0 ? print(" \\hline \n") : print(" \n")                   #new outer multi row

        if( (i-1) % nmrows1 != 0 )
#            (i-1) % nmrows2 == 0 ? print("\\cline{2-14} \n")  : print(" \n")               #new inner multi row
        else
            #nothing
        end

        print("\t\t\t")

        #columns
        for j=1:ncol( df )
            if j == 1
                (i-1) % nmrows1 == 0 ? print("\\multirow{$nmrows1}{*}{", df.ifile[i], "}") : print("")     #new outer multi row
        #    elseif j == 2
        #        (i-1) % nmrows2 == 0 ? print("& \\multirow{$nmrows2}{*}{", df.kF[i], "}") : print(" & ")        #new inter multi row
            else
                if typeof(df[i,j]) == String
                    print(" & ", df[i,j] )
                else
                    if Int(df[i,j]) >= 1000
                        convertedNumber = insertSpaceAfterThreeDigits( Int(df[i,j]) )
                        print(" & ", convertedNumber )
                    else
                        print(" & ", Int(df[i,j]) )
                    end

                end
            end
        end
        (i-1) % nmrows1 == 2 ? print(" \\hspace{-2.5mm} \\Bstrut") : print("")
        (i-1) % nmrows1 == 0 ? print(" \\hspace{-2.5mm} \\Tstrut") : print("")
        print(" \\\\")
    end

    print("\n")
    println( "\t\t\\end{tabular}")
    println( "\t\\label{tab:mylabel}")
    println( "\t\\end{adjustbox}")
    println( "\\end{table}")
    println("\n")
end
