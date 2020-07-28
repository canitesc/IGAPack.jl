using Plots
include("SplineModel.jl")

function nrbplot(nurbs::NURBS, subd::Array)#; light="off", cmap="summer")
    # light = lowercase(light)
    # if !isa(light, String)
    #     error("Parameter for light must be a string (on or off)")
    # end
    #
    # if !isa(cmap, String)
    #     error("Parameter for cmap must be a string")
    # end
    # if (light!="off") || (light=="on")
    #     shade = true
    # else
    #     shade = false
    # end
    subd .+=  1
    #plot the curve or surface
    if length(nurbs.knots)==2
        knt = nurbs.knots
        order = nurbs.order
        evalPtsU = LinRange(knt[1][order[1]],knt[1][end-order[1]+1], subd[1])
        evalPtsV = LinRange(knt[2][order[2]],knt[2][end-order[2]+1], subd[2])
        p = nrbeval(nurbs, [collect(evalPtsU), collect(evalPtsV)])
        plt = Plots.surface!(p[1,:,:], p[2,:,:], p[3,:,:], fc = :viridis)
    elseif length(nurbs.knots)==3
        bnd = nrbextract(nurbs)
        nrbplot(bnd[1], subd[2:3])
        nrbplot(bnd[2], subd[2:3])
        nrbplot(bnd[3], subd[[1,3]])
        nrbplot(bnd[4], subd[[1,3]])
        nrbplot(bnd[5], subd[1:2])
        plt = nrbplot(bnd[6], subd[1:2])
    elseif length(nurbs.knots)==1
        order = nurbs.order
        knt = nurbs.knots
        evalPts = LinRange(knt[1][order[1]],knt[1][end-order[1]+1], subd[1])
        p = nrbeval(nurbs, [collect(evalPts)])
        if any(nurbs.coefs[3,:].!=0)
            plt = Plots.plot!(p[1,:],p[2,:],p[3,:],lw=2, aspect_ratio=:equal)
        else
            plt = Plots.plot!(p[1,:],p[2,:],lw=2, aspect_ratio=:equal)
        end
        return plt
    end
end

function nrbkntplot(nurbs::NURBS, numDim::String="auto")
    nsub = 10

    if length(nurbs.knots)==2 #plot a NURBS surface
        plt = nrbplot(nurbs, [nsub, nsub])

        #plot the knots
        knt1 = unique(nurbs.knots[1][nurbs.order[1]:end-nurbs.order[1]+1])
        knt2 = unique(nurbs.knots[2][nurbs.order[2]:end-nurbs.order[2]+1])
        p1 = nrbeval(nurbs, [knt1, collect(LinRange(knt2[1], knt2[end], nsub))])
        p2 = nrbeval(nurbs, [collect(LinRange(knt1[1], knt1[end], nsub)), knt2])
        if numDim=="auto" && any(nurbs.coefs[3,:,:].!=0)
            numDim="3"
        elseif numDim=="auto" && !any(nurbs.coefs[3,:,:].!=0)
            numDim="2"
        end
        if numDim=="3"
            #surface in a 3D space
            for ii = 1:length(knt1)
                plt = Plots.plot!(p1[1,ii,:], p1[2,ii,:], p1[3,ii,:], linewidth = 2,
                         color=:black, legend=false)
            end
            for ii = 1:length(knt2)
                plt = Plots.plot!(p2[1,:,ii], p2[2,:,ii], p2[3,:,ii], linewidth = 2,
                        color=:black, legend=false)
            end
        else
            #plane surface
            for ii = 1:length(knt1)
                plt = Plots.plot!(p1[1,ii,:], p1[2,ii,:], color=:black, legend=false)
            end
            for ii = 1:length(knt2)
                plt = Plots.plot!(p2[1,ii,:], p2[2,ii,:], color=:black, legend=false)
            end
        end
    elseif length(nurbs.knots)==3 #plot a NURBS volume
        if numDim=="auto"
            numDim="3"
        end
        # Plot the boundaries
        bnd = nrbextract(nurbs)
        for iface = 1:6
            plt = nrbkntplot(bnd[iface], numDim)
        end
    else #plot a NURBS curve
        nrbplot(nurbs, [nsub])
        if numDim=="auto" && any(nurbs.coefs[3,:].!=0)
            numDim="3"
        elseif numDim=="auto" && !any(nurbs.coefs[3,:].!=0)
            numDim="2"
        end
        # And plot the knots
        order = nurbs.order[1]
        p = nrbeval(nurbs, unique(nurbs.knots[1][order:end-order+1]))
        if numDim=="3"
            plt = Plots.scatter!(p[1,:], p[2,:], p[3,:], color=:red, markersize=1,
             shape=:xcross, legend=false)
        else
            plt = Plots.scatter!(p[1,:], p[2,:], color=:red, shape=:xcross,
                    legend=false)
        end
    end
    return plt
end

function nrbctrlplot(nurbs::NURBS, numDim::String="auto")
    if length(nurbs.knots)==3
        nsub = 100
        plt = nrbplot(nurbs, [nsub, nsub, nsub])
        # Plot the control points
        coefs = broadcast(/, nurbs.coefs[1:3,:,:,:], nurbs.coefs[[4],:,:,:])
        coefs = reshape(coefs, 3, :)
        @show coefs
        readline(stdin)
        plt = Plots.scatter!(coefs[1,:], coefs[2,:], coefs[3,:], color=:red,
                markersize=2, legend = false)
        # Plot the control net
        for ii = 1:size(nurbs.coefs, 2)
            for jj = 1:size(nurbs.coefs, 3)
                coefs = reshape(nurbs.coefs[1:3,ii,jj,:], 3, :)
                weights = reshape(nurbs.coefs[[4],ii,jj,:], :)
                xval = coefs[1,:]./weights
                yval = coefs[2,:]./weights
                zval = coefs[3,:]./weights
                plt = Plots.plot!(xval, yval, zval, color = :black,
                    linewidth = 2, linestyle = :dash)
            end
            for kk = 1:size(nurbs.coefs, 4)
                coefs = reshape(nurbs.coefs[1:3,ii,:,kk], 3, :)
                weights = reshape(nurbs.coefs[4,ii,:,kk], :)
                xval = coefs[1,:]./weights
                yval = coefs[2,:]./weights
                zval = coefs[3,:]./weights
                plt = Plots.plot!(xval, yval, zval, color = :black,
                    linewidth = 2, linestyle = :dash)
            end
        end
        for jj = 1:size(nurbs.coefs, 3)
            for kk = 1:size(nurbs.coefs, 4)
                coefs = reshape(nurbs.coefs[1:3, :, jj, kk], 3, :)
                weights = reshape(nurbs.coefs[4, :, jj, kk], :)
                xval = coefs[1,:]./weights
                yval = coefs[2,:]./weights
                zval = coefs[3,:]./weights
                plt = Plots.plot!(xval, yval, zval, color = :black,
                    linewidth = 2, linestyle = :dash)
            end
        end
    elseif length(nurbs.knots)==2
        # Plot a NURBS surface
        nsub = 100
        nrbplot(nurbs, [nsub, nsub])

        # And plot the control net
        for ii = 1:size(nurbs.coefs, 2)
            coefs = reshape(nurbs.coefs[1:3, ii, :], 3, :)
            weights = reshape(nurbs.coefs[4, ii, :], :)
            xval = coefs[1,:]./weights
            yval = coefs[2,:]./weights
            zval = coefs[3,:]./weights
            plt = Plots.plot!(xval, yval, zval, color = :black,
                linewidth = 2, linestyle = :dash)
            plt = Plots.scatter!(xval, yval, zval, color=:red,
                    markersize=2, legend = false)
        end
        for jj = 1:size(nurbs.coefs, 3)
            coefs = reshape(nurbs.coefs[1:3, :, jj], 3, :)
            weights = reshape(nurbs.coefs[4, :, jj], :)
            xval = coefs[1,:]./weights
            yval = coefs[2,:]./weights
            zval = coefs[3,:]./weights
            plt = Plots.plot!(xval, yval, zval, color = :black,
                linewidth = 2, linestyle = :dash)
            plt = Plots.scatter!(xval, yval, zval, color=:red,
                    markersize=2, legend = false)
        end
    else
        # Plot a NURBS curve
        nsub = 1000
        nrbplot(nurbs, [nsub])
        if numDim=="auto" && any(nurbs.coefs[3,:].!=0)
            numDim="3"
        elseif numDim=="auto" && !any(nurbs.coefs[3,:].!=0)
            numDim="2"
        end
        # And plot the control polygon
        coefs = nurbs.coefs[1:3, :]
        weights = nurbs.coefs[4,:]
        xval = coefs[1,:]./weights
        yval = coefs[2,:]./weights
        zval = coefs[3,:]./weights
        if numDim=="3"
            plt = Plots.plot!(xval, yval, zval, color = :black,
                linewidth = 2, linestyle = :dash)
            plt = Plots.scatter!(xval, yval, zval, color=:red, markersize=2,
                legend=false)
        else
            plt = Plots.plot!(xval, yval, color = :black,
                    linewidth = 2, linestyle = :dash)
            plt = Plots.scatter!(xval, yval, color=:red, markersize=5,
                    legend=false)
        end
    end
    return plt
end
