using LinearAlgebra
function bernsteinBasis(u, deg)
    """
    Function returns the values of  Bernstein basis of degree deg at points u
    Algorithm A1.3 in Piegl & Tiller
    xi is a 1D array
    """
    lenU = length(u)
    B = zeros(lenU, deg+1)
    B[:, 1] = ones(lenU)
    u1 = ones(lenU) - u
    u2 = ones(lenU) + u

    for j in 1:deg
        saved = zeros(lenU)
        for k in 1:j
            temp = B[:,k]
            B[:,k] = saved + u1.*temp
            saved  = u2.*temp
        end
        B[:,j+1] = saved
    end
    B = B./(2^deg)

    #calculate the 1st derivative of Bernstein polynomials
    dB = zeros(lenU, deg)
    dB[:,1] = ones(lenU)
    for j in 1:deg-1
        saved = zeros(lenU)
        for k in 0:j-1
            temp = dB[:,k+1]
            dB[:,k+1] = saved + u1.*temp
            saved = u2.*temp
        end
        dB[:,j+1] = saved
    end
    dB = dB./(2^deg)
    dB = hcat(zeros(lenU,1), dB, zeros(lenU, 1))
    dB = (dB[:,1:end-1]-dB[:,2:end])*deg

    #calculate the 2nd derivative of Bernstein polynomials
    if deg>1
        ddB = zeros(lenU, deg-1)
        ddB[:,1] = ones(lenU)
        for j=1:deg-2
            saved = zeros(lenU)
            for k in 0:j-1
                temp = ddB[:,k+1]
                ddB[:,k+1] = saved + u1.*temp
                saved = u2.*temp
            end
            ddB[:,j+1] = saved
        end
        ddB = ddB./(2^deg)
        ddB = hcat(zeros(lenU,2), ddB, zeros(lenU,2))
        ddB = (ddB[:,1:end-2]-2*ddB[:,2:end-1]+ddB[:,3:end])*deg*(deg-1)
    else
        ddB = zeros(lenU,deg+1)
    end
    return B, dB, ddB


end

function bezierExtraction(knot, deg)
    """
    Bezier extraction
    Based on Algroithm 1, from Borden - Isogeometric finite element data
    structures based on Bezier extraction
    """
    m = length(knot)-deg-1;
    a = deg + 1;
    b = a + 1;
    nb = 1;
    nb_final = length(unique(knot))-1
    C = zeros(deg+1,deg+1,nb_final)
    C[:,:,1] = Matrix{Float64}(I, deg+1, deg+1)

    while b<=m
        C[:,:,nb+1] = Matrix{Float64}(I, deg+1, deg+1)
        i = b
        while (b<=m) && (knot[b+1] == knot[b])
            b = b + 1
        end
        multiplicity = b - i + 1
        alphas = zeros(deg-multiplicity)
        if multiplicity < deg
            numerator = knot[b] - knot[a]
            for j in deg:-1:multiplicity+1
                alphas[j-multiplicity] = numerator/(knot[a+j]-knot[a])
            end
            r = deg - multiplicity
            for j in 1:r
                save = r - j + 1
                s = multiplicity + j
                for k in deg+1:-1:s+1
                    alpha = alphas[k-s]
                    C[:,k,nb] = alpha*C[:,k,nb]+(1-alpha)*C[:,k-1,nb]
                end
                if b<=m
                    C[save:save+j,save,nb+1] = C[deg-j+1:deg+1,deg+1,nb]
                end
            end
            nb = nb + 1
            if b <= m
                a = b
                b = b + 1
            end
        elseif multiplicity==deg
            if b <= m
                nb = nb + 1
                a = b
                b = b + 1
            end
        end
    end
    @assert nb_final == nb
    return C, nb
end

"""
Create the extended knot vector (Subsection 4.3.2 in Scott - Isogeometric
data structures based on the Bézier extraction of T-Splines)
"""
function formExtendedKnot(localKnot, p)
    # Repeat the first knot (if needed) so that it appears p+1 times
    firstKnot = localKnot[1]
    indexFirst = findall(localKnot.==firstKnot)
    numRep = length(indexFirst)
    numNewRepFirst = p+1-numRep

    #repeat the last knot (if needed) so that it appears p+1 times
    lastKnot = localKnot[end]
    indexLast = findall(localKnot.==lastKnot)
    numRep = length(indexLast)
    numNewRepLast = p+1-numRep

    #form the extended knot vector
    extendedKnot = vcat(firstKnot*ones(numNewRepFirst), localKnot, lastKnot*ones(numNewRepLast))
    indexFun = numNewRepFirst + 1
    return extendedKnot, indexFun
end

"""
Compute the Bézier extraction operator corresponding to the basis functions of
degree p with local knot vector localKnot
"""
function bezierExtractionLocal(localKnot, p)
    extendedKnot, indexFun = formExtendedKnot(localKnot, p)
    #perform Bézier extraction and return the basis with index numNewRepFirst+1
    C_temp, nb = bezierExtraction(extendedKnot, p)
    IEN, _ = makeIEN(extendedKnot, nb, p)
    C = zeros(p+1, nb)
    for indexSpan = 1:nb
        C_index = findfirst(IEN[indexSpan,:].==indexFun)
        C[:, indexSpan] = C_temp[C_index,:,indexSpan]
    end
    return C
end
