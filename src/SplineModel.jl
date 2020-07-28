include("VecUtils.jl")

struct NURBS
    number::Array{Int64,1}
    coefs::Array{Float64}
    knots::Array{Array{Float64,1}}
    order::Array{Int64,1}
end

function nrbmak(coefs, knots)
    np = [i for i in size(coefs)]
    dim = np[1]

    if length(knots)==3
        if length(np)==3
            np = push!(np,1)
        elseif length(np)==2
            np = append!(np,[1,1])
        end
        #constructing a volume
        number = np[2:4]
        if dim<4
            temp_coefs = coefs
            coefs = repeat([0., 0., 0., 1.], 1, np[2], np[3], np[4])
            coefs[1:dim,:,:,:] = temp_coefs
        end
        uorder = length(knots[1])-np[2]
        vorder = length(knots[2])-np[3]
        worder = length(knots[3])-np[4]
        uknots = sort(knots[1])
        vknots = sort(knots[2])
        wknots = sort(knots[3])
        knots = [uknots, vknots, wknots]
        order = [uorder, vorder, worder]
    elseif length(knots)==2
        if length(np)==2
            np = push!(np,1)
        end
        #constructing a surface
        number = np[2:3]
        if dim<4
            temp_coefs = coefs
            coefs = repeat([0., 0., 0., 1.], 1, np[2], np[3])
            coefs[1:dim,:,:,:] = temp_coefs
        end
        uorder = length(knots[1])-np[2]
        vorder = length(knots[2])-np[3]
        uknots = sort(knots[1])
        vknots = sort(knots[2])
        knots = [uknots, vknots]
        order = [uorder, vorder]
    else
        #constructing a curve
        if length(np)==1
            number = [np[1]]
        else
            number = [np[2]]
        end
        if dim<4
            temp_coefs = coefs
            coefs = repeat([0., 0., 0., 1.], 1, number[1])
            coefs[1:dim,:,:,:] = temp_coefs
        end
        order = [length(knots[1])-number[1]]
        uknots = sort(knots[1])
        knots = [uknots]
    end

    nrb = NURBS(number, coefs, knots, order)
    return nrb
end

function findspan(n::Int64, u, knot::Array{Float64,1})
    if (minimum(u)<knot[1]) || (maximum(u)>knot[end])
        error("Some value is outside the knot span")
    end
    if isa(u, Number)
        s = [0]
    else
        s = similar(u, Int64)
    end
    for j=1:length(u)
        if u[j]==knot[n+2]
            s[j]=n
            continue
        end
        s[j]=findlast(knot.<=u[j])-1
    end
    return s
end

function basisfun(iv::Array{Int64}, uv, p::Integer, U)
    B = zeros(length(uv), p+1)
    N = zeros(p+1)
    for jj = 1:length(uv)
        i = iv[jj]+1 #findspan uses 0-based numbering
        u = uv[jj]
        left = zeros(p+1)
        right = zeros(p+1)
        N[1] = 1

        for j=1:p
            left[j+1] = u - U[i+1-j]
            right[j+1] = U[i+j] - u;
            saved = 0
            for r = 0:j-1
                temp = N[r+1]/(right[r+2]+left[j-r+1])
                N[r+1] = saved + right[r+2]*temp
                saved = left[j-r+1]*temp
            end
            N[j+1] = saved
        end
        B[jj,:] = N
    end
    return B
end

function bspeval(d, c, k, u)
    nu = length(u)
    mc, nc = size(c)
    #readline(stdin)
    s = findspan(nc-1, u, k)
    N = basisfun(s, u, d, k)
    tmp1 = s .- d .+ 1
    p = zeros(mc, nu)
    for i=0:d
        p = p + repeat(N[:,i+1]',mc,1).*c[:,tmp1.+i]
    end
    return p
end

function nrbeval(nurbs::NURBS, tt::AbstractArray)
    foption = true
    degree = nurbs.order .- 1
    if length(nurbs.knots)==3
        #NURBS structure represents a volume
        num1 = nurbs.number[1]
        num2 = nurbs.number[2]
        num3 = nurbs.number[3]

        if isa(tt[1], Array)
            nt1 = length(tt[1])
            nt2 = length(tt[2])
            nt3 = length(tt[3])

            #Evaluate along the w direction
            val = reshape(nurbs.coefs, 4*num1*num2, num3)
            val = bspeval(degree[3], val, nurbs.knots[3], tt[3])
            val = reshape(val, (4, num1, num2, nt3))

            #Evaluate along the v direction
            val = permutedims(val, [1, 2, 4, 3])
            val = reshape(val, 4*num1*nt3, num2)
            val = bspeval(degree[2], val, nurbs.knots[2], tt[2])
            val = reshape(val, 4, num1, nt3, nt2)
            val = permutedims(val, [1, 2, 4, 3])

            #Evaluate along the u direction
            val = permutedims(val, [1, 3, 4, 2])
            val = reshape(val, 4*nt2*nt3, num1)
            val = bspeval(degree[1], val, nurbs.knots[1], tt[1])
            val = reshape(val, 4, nt2, nt3, nt1)
            val = permutedims(val, [1, 4, 2, 3])
            pnts = val
            p = pnts[1:3,:,:,:]
            w = pnts[4,:,:,:,:]'
            if foption
                p = p./repeat(w, 3, 1, 1, 1)
            end
        else
            #Evaluate at scattered points
            #tt[1,:] represents the u direction
            #tt[2,:] represents the v direction
            #tt[3,:] represents the w direction
            st = size(tt)
            if (st[1]!=3 && st[2]==3 && length(st)==2)
                tt = tt'
                st = size(tt)
            end
            nt = prod(st[2:end])
            tt = reshape(tt, (3, nt))
            # Evaluate along the w direction
            val = reshape(nurbs.coefs, 4*num1*num2, num3)
            val = bspeval(degree[3], val, nurbs.knots[3], tt[3,:])
            val = reshape(val, (4, num1, num2, nt))
            # Evaluate along the v direction
            val2 = zeros(4*num1, nt)
            for v=1:nt
                coefs = reshape(val[:,:,:,v], (4*num1, num2))
                val2[:,v] = bspeval(degree[2], coefs, nurbs.knots[2], tt[2,v])
            end
            val2 = reshape(val2, (4, num1, nt))
            # Evaluate along the u direction
            pnts = zeros(4, nt)
            for v = 1:nt
                coefs = reshape(val2[:,:,v], (4, num1))
                pnts[:,v] = bspeval(degree[1], coefs, nurbs.knots[1], tt[1,v])
            end
            w = pnts[4,:]'
            p = pnts[1:3,:]
            if (foption)
                p = p./repeat(w,3,1)
            end

            if length(st)!=2
                w = reshape(w, st[2:end])
                p = reshape(p, tuple((3,)...,st[2:end]...))
            end
        end
    elseif length(nurbs.knots)==2
        # NURBS structure represents a surface
        num1 = nurbs.number[1]
        num2 = nurbs.number[2]
        if isa(tt[1],Array)
            #Evaluate over a [u,v] grid
            # tt[1] represents the u direction
            # tt[2] represents the v direction
            nt1 = length(tt[1])
            nt2 = length(tt[2])

            #Evaluate along the v direction
            val = reshape(nurbs.coefs, 4*num1, num2)
            val = bspeval(degree[2], val, nurbs.knots[2], tt[2])
            val = reshape(val, 4, num1, nt2)

            #Evaluate along the u direction
            val = permutedims(val, [1, 3, 2])
            val = reshape(val, 4*nt2, num1)
            val = bspeval(degree[1], val, nurbs.knots[1], tt[1])
            val = reshape(val, 4, nt2, nt1)
            val = permutedims(val, [1, 3, 2])

            w = val[[4],:,:]
            p = val[1:3, :, :]
            if (foption)
                p = p./repeat(w, 3, 1, 1)
            end
        else
            # Evaluate at scattered points
            # tt[1,:] represents the u direction
            # tt[2,:] represents the v direction
            st = size(tt)
            if (st[1]!=2 && st[2]==2 && length(st)==2)
                tt = tt'
                st = size(tt)
            end
            nt = prod(st[2:end])
            tt = reshape(tt, (3, nt))
            # Evaluate along the v direction
            val = reshape(nurbs.coefs, 4*num1, num2)
            val = bspeval(degree[2], val, nurbs.knots[2], tt[2,:])
            val = reshape(val, (4, num1, nt))
            #Evaluate along the u direction
            pnts = zeros(4, nt)
            for v = 1:nt
                coefs = reshape(val[:,:,v], (4, num1))
                pnts[:,v] = bspeval(degree[1], coefs, nurbs.knots[1], tt[1,v])
            end
            w = pnts[4,:]'
            p = pnts[1:3,:]
            if (foption)
                p = p./repeat(w,3,1)
            end
            if (length(st)!=2)
                w = reshape(w, st[2:end])
                p = reshape(p, tuple((3,)...,st[2:end]...))
            end
        end
    else
        if isa(tt[1], Array) && length(tt)==1
            tt =tt[1]
        end
        st = size(tt)
        val = bspeval(nurbs.order[1]-1, nurbs.coefs, nurbs.knots[1], [tt...])
        w = val[4,:]'
        p = val[1:3,:]
        if foption
            p = p./repeat(w,3,1)
        end
        if (st[1]!=1 || length(st)!=2)
            w = reshape(w, st)
            p = reshape(p, tuple((3,)...,st...))
        end
    end
    return p
end

function nrbextract(srf::NURBS)
    for idim = 1:length(srf.knots)
        ord = srf.order[idim]
        if (srf.knots[idim][1]!=srf.knots[idim][ord] ||
            srf.knots[idim][end]!=srf.knots[idim][end-ord+1])
            error("nrbextract: only working with open knot vectors")
        end
    end

    if (length(srf.knots)==2)
        crvs = Array{NURBS}(undef, 4)
        for ind=1:2
            ind2 = mod(ind,2)+1
            bnd1 = (ind-1)*2+1
            bnd2 = (ind-1)*2+2
            if ind==1
                coefs1 = srf.coefs[:,1,:]
                coefs2 = srf.coefs[:,end,:]
            elseif ind==2
                coefs1 = srf.coefs[:,:,1]
                coefs2 = srf.coefs[:,:,end]
            end
            crvs[bnd1] = nrbmak(coefs1, srf.knots[ind2])
            crvs[bnd2] = nrbmak(coefs2, srf.knots[ind2])
        end
    elseif (length(srf.knots)==3)
        crvs = Array{NURBS}(undef, 6)
        for ind=1:3
            inds = setdiff(1:3, ind)
            bnd1 = (ind - 1)*2 + 1
            bnd2 = (ind - 1)*2 + 2
            if ind==1
                coefs1 = srf.coefs[:,1,:,:]
                coefs2 = srf.coefs[:,end,:,:]
            elseif ind==2
                coefs1 = srf.coefs[:,:,1,:]
                coefs2 = srf.coefs[:,:,end,:]
            elseif ind==3
                coefs1 = srf.coefs[:,:,:,1]
                coefs2 = srf.coefs[:,:,:,end]
            end
            crvs[bnd1] = nrbmak(coefs1, [srf.knots[inds[1]], srf.knots[inds[2]]])
            crvs[bnd2] = nrbmak(coefs2, [srf.knots[inds[1]], srf.knots[inds[2]]])
        end
    else
        error("The entity is not a surface nor a volume")
    end
    return crvs
end

function nrbline(p1, p2)
    coefs = [zeros(3,2); ones(1,2)]
    coefs[1:length(p1),1] .= p1
    coefs[1:length(p2),2] .= p2
    line = nrbmak(coefs, [[0,0,1,1]])
    return line
end

function nrbline()
    p1 = [0,0]
    p2 = [1,0]
    line = nrbline(p1, p2)
    return line
end

function nrbrect(width, height)
    coefs =[0. width width width width 0. 0. 0.;
          0. 0. 0. height height height height 0.;
          0. 0. 0. 0. 0. 0. 0. 0.;
          1. 1. 1. 1. 1. 1. 1. 1.]
    knots = [0, 0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1, 1]
    curve = nrbmak(coefs, [knots])
    return curve
end

function nrbrect(size)
    return nrbrect(size, size)
end

function nrbrect()
    return nrbrect(1.)
end

function nrbcirc(radius::Real, center::Vector, sang::Real, eang::Real)
    sweep = eang - sang
    if sweep < 0
        sweep = 2*pi + sweep
    end

    if abs(sweep) <= pi/2
        narcs = 1   # number of arc segments
        knots = [0, 0, 0, 1, 1, 1]
    elseif abs(sweep) <= pi
        narcs = 2
        knots = [0, 0, 0, 0.5, 0.5, 1, 1, 1]
    elseif abs(sweep) <= 3*pi/2
        narcs = 3
        knots = [0, 0, 0, 1/3, 1/3, 2/3, 2/3, 1, 1, 1]
    else
        narcs = 4
        knots = [0, 0, 0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1, 1, 1]
    end

    dsweep = sweep/(2*narcs)  #arc segment sweep angle/2

    # determine middle control point and weight
    wm = cos(dsweep)
    x  = radius*wm
    y  = radius*sin(dsweep)
    xm = x+y*tan(dsweep)

    # arc segment control points
    ctrlpt = [ x wm*xm x;    # w*x - coordinate
              -y 0     y;    # w*y - coordinate
               0 0     0;    # w*z - coordinate
               1 wm    1]    # w   - coordinate

    # build up complete arc from rotated segments
    coefs = zeros(4,2*narcs+1)   # nurbs control points of arc
    xx = vecrotz(sang + dsweep)
    coefs[:,1:3] = xx*ctrlpt;     # rotate to start angle
    xx = vecrotz(2*dsweep);
    for n = 2:narcs
       m = 2*n.+[0, 1]
       coefs[:,m] = xx*coefs[:,m.-2]
    end

    # vectrans arc if necessary
    if !isempty(center)
      xx = vectrans(center)
      coefs = xx*coefs
    end

    curve = nrbmak(coefs,[knots])
    return curve
end

function nrbcirc(radius::Real, center::Vector)
    sang = 0.
    eang = 2*pi
    curve = nrbcirc(radius, center, sang, eang)
    return curve
end

function nrbcirc(radius::Real)
    center = []
    curve = nrbcirc(radius, center)
    return curve
end

function nrbcirc()
    radius = 1
    curve = nrbcirc(radius)
    return curve
end

function bspkntins(d, c, k, u)
    tolEq = 1e-10
    mc, nc = size(c)
    sort!(u)
    nu = length(u)
    nk = length(k)

    ic = zeros(mc, nc+nu)
    ik = zeros(nk+nu)

    n = nc - 1
    r = nu - 1

    m = n + d + 1
    a = findspan(n, u[1], k)[1]
    b = findspan(n, u[r+1], k)[1]
    b += 1

    ic[:, 1:a-d+1] = c[:, 1:a-d+1]
    ic[:, b+nu:nc+nu] = c[:, b:nc]

    ik[1:a+1] = k[1:a+1]
    ik[b+d+nu+1:m+nu+1] = k[b+d+1:m+1]

    ii = b + d - 1
    ss = ii + nu

    for jj = r:-1:0
        ind = (a+1):ii
        ind = ind[u[jj+1].<=k[ind.+1]]
        ic[:, ind.+ss.-ii.-d] = c[:,ind.-d]
        ik[ind.+ss.-ii.+1] = k[ind.+1]
        ii = ii - length(ind)
        ss = ss - length(ind)

        ic[:,ss-d] = ic[:, ss-d+1]
        for l=1:d
            ind = ss - d + l
            alfa = ik[ss+l+1] - u[jj+1]
            if abs(alfa) < tolEq
                ic[:,ind] = ic[:, ind+1]
            else
                alfa = alfa/(ik[ss+l+1]-k[ii-d+l+1])
                tmp = (1-alfa) * ic[:, ind+1]
                ic[:,ind] = alfa*ic[:,ind] + tmp
            end
        end
        ik[ss+1] = u[jj+1]
        ss = ss-1
    end
    return ic, ik
end

function nrbkntins(nurbs::NURBS, iknots::AbstractArray)
    degree = nurbs.order .- 1
    fmax(x,y) = any(y.>maximum(x))
    fmin(x,y) = any(y.<minimum(x))
    for i=1:length(nurbs.knots)
        if any(fmax(nurbs.knots[i], iknots[i])) || any(fmin(nurbs.knots[i], iknots[i]))
            error("Trying to insert a knot outside the interval of definition")
        end
    end
    knots = Array{Array}(undef, length(nurbs.knots))
    if length(nurbs.knots)==3
        # NURBS represents a volume
        num1 = nurbs.number[1]
        num2 = nurbs.number[2]
        num3 = nurbs.number[3]

        # Insert knots along the w direction
        if isempty(iknots[3])
            coefs = nurbs.coefs
            knots[3] = nurbs.knots[3]
        else
            coefs = reshape(nurbs.coefs, 4*num1*num2, num3)
            coefs, knots[3] = bspkntins(degree[3], coefs, nurbs.knots[3], iknots[3])
            num3 = size(coefs, 2)
            coefs = reshape(coefs, 4, num1, num2, num3)
        end

        # Insert knots along the v direction
        if isempty(iknots[2])
            knots[2] = nurbs.knots[2]
        else
            coefs = permutedims(coefs, [1, 2, 4, 3])
            coefs = reshape(coefs, 4*num1*num3, num2)
            coefs, knots[2] = bspkntins(degree[2], coefs, nurbs.knots[2], iknots[2])
            num2 = size(coefs, 2)
            coefs = reshape(coefs, 4, num1, num3, num2)
            coefs = permutedims(coefs, [1, 2, 4, 3])
        end

        #Insert knots along the u direction
        if isempty(iknots[1])
            knots[1] = nurbs.knots[1]
        else
            coefs = permutedims(coefs, [1, 3, 4, 2])
            coefs = reshape(coefs, 4*num2*num3, num1)
            coefs, knots[1] = bspkntins(degree[1], coefs, nurbs.knots[1], iknots[1])
            coefs = reshape(coefs, 4, num2, num3, size(coefs, 2))
            coefs = permutedims(coefs, [1, 4, 2, 3])
        end
    elseif length(nurbs.knots)==2
        # NURBS rrepresents a surface
        num1 = nurbs.number[1]
        num2 = nurbs.number[2]

        # Insert knots along the v direction
        if isempty(iknots[2])
            coefs = nurbs.coefs
            knots[2] = nurbs.knots[2]
        else
            coefs = reshape(nurbs.coefs, 4*num1, num2)
            coefs, knots[2] = bspkntins(degree[2], coefs, nurbs.knots[2], iknots[2])
            num2 = size(coefs, 2)
            coefs = reshape(coefs, 4, num1, num2)
        end

        # Insert knots along the u directions
        if isempty(iknots[1])
            knots[1] = nurbs.knots[1]
        else
            coefs = permutedims(coefs, [1, 3, 2])
            coefs = reshape(coefs, 4*num2, num1)
            coefs, knots[1] = bspkntins(degree[1], coefs, nurbs.knots[1], iknots[1])
            coefs = reshape(coefs, 4, num2, size(coefs,2))
            coefs = permutedims(coefs, [1, 3, 2])
        end
    else
        # NURBS represents a curve
        if isempty(iknots[1])
            coefs = nurbs.coefs
            knots = nurbs.knots
        else
            coefs, knots[1] = bspkntins(degree[1], nurbs.coefs, nurbs.knots[1], iknots[1])
        end
    end

    #construct the new NURBS
    inurbs = nrbmak(coefs, knots)
    return inurbs
end

function bspdegelev(d, c, k, t)
    mc, nc = size(c)
    ic = zeros(mc, nc*(t+1))
    n = nc - 1
    bezalfs = zeros(d+1, d+t+1)
    bpts = zeros(mc, d+1)
    ebpts = zeros(mc, d+t+1)
    Nextbpts = zeros(mc, d+1)
    alfs = zeros(d)

    m = n + d + 1
    ph = d + t
    ph2 = floor(Int, ph/2)
    # Compute Bezier degree elevation coefficients
    bezalfs[1,1] = 1
    bezalfs[d+1, ph+1] = 1

    for i=1:ph2
        inv = 1/binomial(ph, i)
        mpi = min(d,i)

        for j=max(0, i-t):mpi
            bezalfs[j+1, i+1] = inv*binomial(d,j)*binomial(t, i-j)
        end
    end

    for i=ph2+1:ph-1
        mpi = min(d,i)
        for j=max(0,i-t):mpi
            bezalfs[j+1, i+1] = bezalfs[d-j+1, ph-i+1]
        end
    end

    mh = ph
    kind = ph+1
    r = -1
    a = d
    b = d + 1
    cind = 1
    ua = k[1]

    ic[1:mc,1] = c[1:mc,1]
    ik = ua.*ones(ph+1)

    # Initialize the first Bezier segment
    #bpts[1:mc, 1:d+1] = c[1:mc, 1:d+1]
    bpts = copy(c)

    # Big loop through knot vector
    while b<m
        i=b
        while b<m && k[b+1] == k[b+2]
            b += 1
        end
        mul = b - i + 1
        mh = mh + mul + t
        ub = k[b+1]
        oldr = r
        r = d - mul

        # Insert knot u[b] r times
        if oldr > 0
            lbz = floor(Int, (oldr+2)/2)
        else
            lbz = 1
        end

        if r>0
            rbz = ph - floor(Int, (r+1)/2)
        else
            rbz = ph
        end

        if r>0
            # Insert knot to get Bezier segment
            numer = ub - ua
            for q = d:-1:mul+1
                alfs[q-mul] = numer / (k[a+q+1]-ua)
            end

            for j = 1:r
                save = r-j
                s = mul + j

                for q = d:-1:s
                    for ii = 0:mc-1
                        tmp1 = alfs[q-s+1]*bpts[ii+1, q+1]
                        tmp2 = (1-alfs[q-s+1])*bpts[ii+1, q]
                        bpts[ii+1, q+1] = tmp1 + tmp2
                    end
                end

                Nextbpts[:, save+1] = bpts[:, d+1]
            end
        end
        # End of insert knot

        #Degree elevate Bezier
        for i=lbz:ph
            ebpts[:,i+1] = zeros(mc)
            mpi = min(d,i)
            for j=max(0,i-t):mpi
                for ii=0:mc-1
                    tmp1 = ebpts[ii+1, i+1]
                    tmp2 = bezalfs[j+1, i+1]*bpts[ii+1, j+1]
                    ebpts[ii+1, i+1] = tmp1 + tmp2
                end
            end
        end
        # End of degree elevating Bezier

        if oldr > 1
            # Must remove knot u=k[a] oldr times
            first = kind - 2
            last = kind
            den = ub - ua
            bet = floor(Int, (ub-ik[kind])/den)

            # Knot removal loop
            for tr = 1:oldr-1
                i = first
                j = last
                kj = j - kind + 1
                while j-i > tr
                    # Loop and compute the new control points for one removal step
                    if i < cind
                        alf = (ub - ik[i+1])/(ua - ik[i+1])
                        tmp1 = alf.*ic[:,i+1]
                        tmp2 = (1-alf).*ic[:, i]
                        ic[:,i+1] = tmp1 + tmp2
                    end
                    if j >= lbz
                        if j-tr <= kind - ph + oldr
                            gam = (ub-ik[j-tr+1])/den
                            tmp1 = gam.*ebpts[:, kj+1]
                            tmp2 = (1-gam).*ebpts[:, kj+2]
                            ebpts[:, kj+1] = tmp1 + tmp2
                        else
                            tmp1 = bet.*ebpts[:, kj+1]
                            tmp2 = (1-bet).*ebpts[:, kj+2]
                            ebpts[:, kj+1] = tmp1 + tmp2
                        end
                    end
                    i += 1
                    j -= 1
                    kj -= 1
                end
                first -= 1
                last += 1
            end
        end
        # End of removing knot n=k[a]

        # Load the knot ua
        if a!=d
            for i=0:ph-oldr-1
                push!(ik, ua)
                #ik[kind+1] = ua
                kind += 1
            end
        end

        for j = lbz:rbz
            for ii = 0:mc-1
                ic[ii+1, cind+1] = ebpts[ii+1, j+1]
            end
            cind += 1
        end

        if b < m
            # Setup for next pass through loop
            bpts[:,1:r] = Nextbpts[:, 1:r]
            bpts[:,r+1:d+1] = c[:, b-d+r+1:b+1]
            a = b
            b += 1
            ua = ub
        else
            for i=0:ph
                #ik[kind+i+1] = ub
                push!(ik, ub)
            end
        end
    end
    #End big while loop
    ic = ic[:, 1:cind]
    return ic, ik
end

function nrbdegelev(nurbs::NURBS, ntimes::Array{Int,1})
    degree = nurbs.order .- 1
    knots = Array{Array}(undef, length(nurbs.knots))
    if length(nurbs.knots)==3
        # NURBS represents a volume
        _, num1, num2, num3 = size(nurbs.coefs)

        # Degree elevate along the w direction
        if ntimes[3] == 0
            coefs = nurbs.coefs
            knots[3] = nurbs.knots[3]
        else
            coefs = reshape(nurbs.coefs, 4*num1*num2,num3)
            coefs, knots[3] = bspdegelev(degree[3], coefs, nurbs.knots[3], ntimes[3])
            num3 = size(coefs, 2)
            coefs = reshape(coefs, 4, num1, num2, num3)
        end

        # Degree elevate along the v direction
        if ntimes[2] == 0
            knots[2] = nurbs.knots[2]
        else
            coefs = permutedims(coefs, [1, 2, 4, 3])
            coefs = reshape(coefs, 4*num1*num3, num2)
            coefs, knots[2] = bspdegelev(degree[2], coefs, nurbs.knots[1], ntimes[1])
            num2 = size(coefs, 2)
            coefs = reshape(coefs, 4, num1, num3, num2)
            coefs = permutedims(coefs, [1, 2, 4, 3])
        end

        # Degree elevate along the u direction
        if ntimes[1] == 0
            knots[1] = nurbs.knots[1]
        else
            coefs = permutedims(coefs, [1, 3, 4, 2])
            coefs = reshape(coefs, 4*num2*num3, num1)
            coefs, knots[1] = bspdegelev(degree[1], coefs, nurbs.knots[1], ntimes[1])
            coefs = reshape(coefs, 4, num2, num3, size(coefs, 2))
            coefs = permutedims(coefs, [1, 4, 2, 3])
        end
    elseif length(nurbs.knots)==2
        _, num1, num2 = size(nurbs.coefs)

        # Degree elevate along the v direction
        if ntimes[2] == 0
            coefs = nurbs.coefs
            knots[2] = nurbs.knots[2]
        else
            coefs = reshape(nurbs.coefs, 4*num1, num2)
            coefs, knots[2] = bspdegelev(degree[2], coefs, nurbs.knots[2], ntimes[2])
            num2 = size(coefs, 2)
            coefs = reshape(coefs, 4, num1, num2)
        end

        # Degree elevate along the u direction
        if ntimes[1] == 0
            knots[1] = nurbs.knots[1]
        else
            coefs = permutedims(coefs, [1, 3, 2])
            coefs = reshape(coefs, 4*num2, num1)
            coefs, knots[1] = bspdegelev(degree[1], coefs, nurbs.knots[1], ntimes[1])
            coefs = reshape(coefs, 4, num2, size(coefs, 2))
            coefs = permutedims(coefs, [1, 3, 2])
        end
    else
        # NURBS represents a curve
        if isempty(ntimes)
            coefs = nurbs.coefs
            knots = nurbs.knots
        else
            coefs, knots[1] = bspdegelev(degree[1], nurbs.coefs, nurbs.knots[1], ntimes[1])
        end
    end

    #construct new NURBS
    inurbs = nrbmak(coefs, knots)
    return inurbs
end

function nrbruled(crv1::NURBS, crv2::NURBS)
    # Ensure both curves have a common degree
    d = max(crv1.order, crv2.order)[1]
    crv1 = nrbdegelev(crv1, [d - crv1.order[1]])
    crv2 = nrbdegelev(crv2, [d - crv2.order[1]])

    # Merge the knot vectors to obtain a common knot vector
    k1 = crv1.knots[1]
    k2 = crv2.knots[1]

    ku = unique(vcat(k1, k2))
    n = length(ku)

    ka = Array{Float64,1}(undef,0)
    kb = Array{Float64,1}(undef,0)

    for i = 1:n
        i1 = length(findall(k1.==ku[i]))
        i2 = length(findall(k2.==ku[i]))
        m = max(i1, i2)
        ka = vcat(ka, ku[i].*ones(m-i1))
        kb = vcat(kb, ku[i].*ones(m-i2))
    end
    crv1 = nrbkntins(crv1, [ka])
    crv2 = nrbkntins(crv2, [kb])

    coefs = cat(crv1.coefs, crv2.coefs, dims=3)
    srf = nrbmak(coefs, [crv1.knots[1], [0,0,1,1]])
    return srf
end

function nrbtform(nurbs::NURBS, tmat)
    if length(nurbs.knots)==2
        # NURBS is a surface
        dim, nu, nv = size(nurbs.coefs)
        coefs = reshape(tmat*reshape(nurbs.coefs, dim, nu*nv), dim, nu, nv)
    elseif length(nurbs.knots)==3
        # NURBS is a volume
        dim, nu, nv, nw = size(nurbs.coefs)
        coefs = reshape(tmat*reshape(nurbs.coefs, dim, nu*nv*nw), dim, nu, nv, nw)
    else
        # NURBS is a curve
        coefs = tmat*nurbs.coefs
    end
    tnurbs = nrbmak(coefs, nurbs.knots)
    return tnurbs
end

function nrbtransp(srf)
    if length(srf.knots)==1
        error("A NURBS curve cannot be transposed")
    elseif length(srf.knots)==3
        error("The transposition of NURBS volumes has not been implemented")
    end
    tsrf = nrbmak(permutedims(srf.coefs, [1,3,2]), reverse(srf.knots))
    return tsrf
end

function nrb4surf(p11, p12, p21, p22)
    coefs = cat(zeros(3,2,2), ones(1,2,2), dims=1)
    coefs[1:length(p11),1,1] = p11[:]
    coefs[1:length(p12),2,1] = p12[:]
    coefs[1:length(p21),1,2] = p21[:]
    coefs[1:length(p22),2,2] = p22[:]

    knots = [[0,0,1,1], [0,0,1,1]]
    srf = nrbmak(coefs, knots)
    return srf
end


function nrbcoons(u1::NURBS, u2::NURBS, v1::NURBS, v2::NURBS)
    tolEq = 1e-10
    if maximum(abs.(nrbeval(u1,u1.knots[1][[1]]) - nrbeval(v1,v1.knots[1][[1]])))>tolEq ||
            maximum(abs.(nrbeval(u1,u1.knots[1][[end]]) - nrbeval(v2,v2.knots[1][[1]])))>tolEq ||
            maximum(abs.(nrbeval(u2,u2.knots[1][[1]]) - nrbeval(v1,v1.knots[1][[end]])))>tolEq ||
            maximum(abs.(nrbeval(u2,u2.knots[1][[end]]) - nrbeval(v2,v2.knots[1][[end]])))>tolEq
        error("The four curves do not define a closed boundary")
    end

    r1 = nrbruled(u1, u2)
    r2 = nrbtransp(nrbruled(v1, v2))
    t = nrb4surf(u1.coefs[:,1], u1.coefs[:,end], u2.coefs[:,1], u2.coefs[:,end])

    # Raise all surfaces to a common degree
    du = max(r1.order[1], r2.order[1], t.order[1])
    dv = max(r1.order[2], r2.order[2], t.order[2])
    r1 = nrbdegelev(r1, [du - r1.order[1], dv - r1.order[2]])
    r2 = nrbdegelev(r2, [du - r2.order[1], dv - r2.order[2]])
    t = nrbdegelev(t, [du - t.order[1], dv - t.order[2]])

    # Merge the knots vectors, to obtain a commmon knot vector

    #Uknots
    k1 = r1.knots[1]
    k2 = r2.knots[1]
    k3 = t.knots[1]
    k = unique(vcat(k1, k2, k3))
    n = length(k)
    kua = Array{Float64,1}(undef,0)
    kub = Array{Float64,1}(undef,0)
    kuc = Array{Float64,1}(undef,0)

    for i=1:n
        i1 = length(findall(k1.==k[i]))
        i2 = length(findall(k2.==k[i]))
        i3 = length(findall(k3.==k[i]))
        m = max(i1, i2, i3)
        kua = vcat(kua, k[i].*ones(m-i1))
        kub = vcat(kub, k[i].*ones(m-i2))
        kuc = vcat(kuc, k[i].*ones(m-i3))
    end

    # Vknots
    k1 = r1.knots[2]
    k2 = r2.knots[2]
    k3 = t.knots[2]

    k = unique(vcat(k1, k2, k3))
    n = length(k)

    kva = Array{Float64,1}(undef,0)
    kvb = Array{Float64,1}(undef,0)
    kvc = Array{Float64,1}(undef,0)

    for i=1:n
        i1 = length(findall(k1.==k[i]))
        i2 = length(findall(k2.==k[i]))
        i3 = length(findall(k3.==k[i]))
        m = max(i1, i2, i3)
        kva = vcat(kva, k[i].*ones(m-i1))
        kvb = vcat(kvb, k[i].*ones(m-i2))
        kvc = vcat(kvc, k[i].*ones(m-i3))
    end

    r1 = nrbkntins(r1, [kua, kva])
    r2 = nrbkntins(r2, [kub, kvb])
    t = nrbkntins(t, [kuc, kvc])

    # combine coefficients to construct the Coons surface
    coefs = zero(r1.coefs)
    coefs[1,:,:] = r1.coefs[1,:,:] + r2.coefs[1,:,:] - t.coefs[1,:,:]
    coefs[2,:,:] = r1.coefs[2,:,:] + r2.coefs[2,:,:] - t.coefs[2,:,:]
    coefs[3,:,:] = r1.coefs[3,:,:] + r2.coefs[3,:,:] - t.coefs[3,:,:]
    coefs[4,:,:] = r1.coefs[4,:,:] + r2.coefs[4,:,:] - t.coefs[4,:,:]
    srf = nrbmak(coefs, r1.knots)
end

function nrbrevolve(curve::NURBS, pnt, vec, theta)
    if length(curve.knots)==3
        error("The function nrbrevolve is not yet ready to create volumes")
    end

    # Translate the center point to the origin
    if isempty(pnt)
        pnt = zeros(3)
    end

    if length(pnt)!=3
        error("All point and vecto coordinates must be 3D")
    end

    # Translate and rotate the original curve or surface into alignment with the z-axis
    T = vectrans(-pnt)
    angx = vecangle(vec[1], vec[3])
    RY = vecroty(-angx[1])
    vectmp = RY*vcat(vecnorm(vec), 1.)
    angy = vecangle(vectmp[2], vectmp[3])
    RX = vecrotx(angy[1])
    curve = nrbtform(curve, RX*RY*T)

    # Construct an arc
    arc = nrbcirc(1.0, [], 0.0, theta)
    if length(curve.knots)==2
        # Construct the revolved volume
        coefs = zeros(4, arc.number[1], curve.number[1], curve.number[2])
        angle = vecangle(curve.coefs[2,:,:], curve.coefs[1,:,:])
        radius = dropdims(vecmag(curve.coefs[1:2, :, :]), dims=1)
        for i=1:curve.number[1]
            for j=1:curve.number[2]
                coefs[:,:,i,j] = vecrotz(angle[i,j])*vectrans([0., 0., curve.coefs[3,i,j]])*
                            vecscale([radius[i,j], radius[i,j]])*arc.coefs
                coefs[4,:,i,j] = coefs[4, :, i, j]*curve.coefs[4,i,j]
            end
        end
        surf = nrbmak(coefs, vcat(arc.knots, curve.knots))
    else
        # Construct the revolved surface
        coefs = zeros(4, arc.number[1], curve.number[1])
        angle = vecangle(curve.coefs[2,:], curve.coefs[1,:])
        radius = vecmag(curve.coefs[1:2, :])
        for i = 1:curve.number[1]
            coefs[:, :, i] = vecrotz(angle[i])*vectrans([0., 0., curve.coefs[3,i]])*
                        vecscale([radius[i], radius[i]])*arc.coefs
            coefs[4, :, i] = coefs[4, :, i]*curve.coefs[4, i]
        end
        surf = nrbmak(coefs, vcat(arc.knots, curve.knots))
    end

    # Rotate and vectrans the surface back into position
    T = vectrans(pnt)
    RX = vecrotx(-angy[1])
    RY = vecroty(angx[1])
    surf = nrbtform(surf, T*RY*RX)
    return surf
end

function nrbrevolve(curve::NURBS, pnt, vec)
    surf = nrbrevolve(curve, pnt, vec, 2*pi)
    return surf
end

#end
function kntrefine(knots, n_sub, degree, regularity)
    if length(n_sub)!=length(degree) || length(n_sub)!=length(regularity) ||
        length(n_sub)!=length(knots)
        error("kntrefine: n_sub, degree and regularity must have the same length
            as the number of knot vector")
    end
    rknots = Array{Array}(undef, length(knots))
    new_knots = Array{Array}(undef, length(knots))
    zeta = Array{Array}(undef, length(knots))
    for idim = 1:length(n_sub)
        min_mult = degree[idim] - regularity[idim]
        z = unique(knots[idim])
        nz = length(z)
        deg = sum(knots[idim].==z[1])-1
        rknots[idim] = z[ones(Int, degree[idim]+1)]
        new_knots[idim] = [];
        for ik = 2:nz
            insk = LinRange(z[ik-1], z[ik], n_sub[idim]+2)
            insk = vec(repeat(insk[2:end-1], min_mult,1))
            old_mult = sum(knots[idim].==z[ik])
            mult = max(min_mult, degree[idim] - deg + old_mult)

            rknots[idim] = vcat(rknots[idim], insk, z[ik*ones(Int, mult)])
            if mult>=old_mult
                new_knots[idim] = vcat(new_knots[idim], insk, z[ik*ones(Int, mult-old_mult)])
            else
                new_knots[idim] = vcat(new_knots[idim], insk)
            end
        end
        zeta[idim] = unique(rknots[idim])
    end
    return rknots, zeta, new_knots
end

function nrbsquare(corner, lengthx, lengthy, degree, nsubdiv)
    if isempty(corner)
        corner = [0,0]
    end

    if length(degree)==1
        deg = [degree, degree]
    elseif length(degree)==2
        deg = degree
    else
        error("The degree vector should provide the degree in each direction (one or two values)")
    end

    if length(nsubdiv)==1
        nsub = [nsubdiv, nsubdiv]
    elseif length(nsubdiv)==2
        nsub = nsubdiv
    else
        error("The nsubdiv vector should provide the number of intervals in each direction (one or two values)")
    end

    srf = nrb4surf(corner, corner+[lengthx, 0], corner+[0, lengthy], corner+[lengthx, lengthy])
    srf = nrbdegelev(srf, deg-[1, 1])
    _, _, new_knots = kntrefine(srf.knots, nsub.-1, deg, deg - [1, 1] )
    srf = nrbkntins(srf, new_knots)
    return srf
end
