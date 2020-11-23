## Definition of projective transformation
"""
Projective transformation.
"""
struct Projective
    H::SArray{Tuple{3,3},Float64,2,9}
end

"""
Definition of projective transformation.
"""
function (p::Projective)(x)
    x_ = SA[x[1],x[2],1.0]
    x′_1, x′_2, x′_3 = p.H*x_
    return SA[x′_1/x′_3, x′_2/x′_3]
end

"""
Inverse of projective transformation.
"""
function Base.inv(p::Projective)
    Projective(inv(p.H))
end

"""
Generate matrix for projective transformation,
which satisfy (0,0)↦p00, (1,0)↦p10, (0,1)↦p01, (1,1)↦p11.
"""
function _projectivematrix(p00,p01,p10,p11)
    k1,k2 = hcat(p11-p10,p11-p01)\(p10+p01-p00-p11)
    v1 = (k1+1)*p10-p00
    v2 = (k2+1)*p01-p00
    h = hcat(v1,v2,p00)
    H = vcat(h,[k1,k2,1]')
    return H
end

"""
Generate projective transformation,
which satisfy (0,0)↦p00, (1,0)↦p10, (0,1)↦p01, (1,1)↦p11.
"""
function Projective(p00,p01,p10,p11)
    H = _projectivematrix(p00,p01,p10,p11)
    p = Projective(H)
    return p
end

"""
Generate projective transformation,
which satisfy p00↦q00, p10↦q10, p01↦q01, p11↦q11.
"""
function Projective(p00,p01,p10,p11,q00,q01,q10,q11)
    Hp = _projectivematrix(p00,p01,p10,p11)
    Hq = _projectivematrix(q00,q01,q10,q11)
    H = Hq*inv(Hp)
    p = Projective(H)
    return p
end

## Tests for projective transformation
# p00 = [0,0]+rand(2)
# p01 = [0,1]+rand(2)
# p10 = [1,0]+rand(2)
# p11 = [1,1]+rand(2)

# q00 = [0,0]+rand(2)
# q01 = [0,1]+rand(2)
# q10 = [1,0]+rand(2)
# q11 = [1,1]+rand(2)

# p = Projective(p00,p01,p10,p11,q00,q01,q10,q11)
# norm(p(p00)-q00)
# norm(p(p01)-q01)
# norm(p(p10)-q10)
# norm(p(p11)-q11)
