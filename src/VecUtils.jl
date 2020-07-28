function vecrotx(angle::Real)
    sn = sin(angle)
    cn = cos(angle)
    rx = [1 0 0 0; 0 cn -sn 0; 0 sn cn 0; 0 0 0 1]
    return rx
end

function vecroty(angle::Real)
    sn = sin(angle)
    cn = cos(angle)
    ry = [cn 0 sn 0; 0 1 0 0; -sn 0 cn 0; 0 0 0 1]
    return ry
end

function vecrotz(angle::Real)
    sn = sin(angle)
    sn = sin(angle)
    cn = cos(angle)
    rz = [ cn -sn 0 0; sn cn 0 0; 0 0 1 0; 0 0 0 1]
    return rz
end

function vecangle(num, den)
    if length(num)==1
        ang = [atan(num, den)]
    else
        ang = zero(num)
        for i=1:length(num)
            ang[i] = atan(num[i], den[i])
        end
    end
    index = findall(ang .< 0.0)
    ang[index] = 2*pi .+ ang[index]
    return ang
end

function vecmag(vec)
    mag  = sqrt.(sum(vec.^2, dims=1))
    return mag
end

function vecnorm(vec)
    nvec = vec./repeat(sqrt.(sum(vec.^2, dims=1)), outer=tuple(vcat(size(vec, 1), ones(Int, ndims(vec)-1))...))
    return nvec
end

function vecscale(vector::Vector)
    s = vcat(vector, [0, 0])
    ss = [s[1] 0 0 0; 0 s[2] 0 0; 0 0 s[3] 0; 0 0 0 1]
    return ss
end

function vectrans(vector::Vector)
    v = vcat(vector,[0, 0])
    dd = [1 0 0 v[1]; 0 1 0 v[2]; 0 0 1 v[3]; 0 0 0 1]
    return dd
end
