using Interpolations
using DifferentialEquations
#import GR
#GR.opengks()
#GR.openws(10, "", 210)
#GR.activatews(10)
import Plots

# Needed Packages
if false
    Pkg.add("Interpolations")
    Pkg.add("DifferentialEquations")
    Pkg.add("Plots")
    Pkg.add("GR")
end


# Versions check
if false
    println("computing versions...")
    versions = Pkg.installed()
    
    assert(VERSION == v"0.6.2")
    assert(versions["Interpolations"] == v"0.7.3")
    assert(versions["DifferentialEquations"] == v"4.1.0")
    assert(versions["Plots"] == v"0.15.0")
    assert(versions["GR"] == v"0.26.0")
end


#=

Dispersion for zero order constitutive approximation
kzZ = acos(1 + Z^2((1 - cos(kxX))/X^2 - (1 - cos(ktT))/T^2))

=#

#=

Nt cubic cells on t direction, periodic condition, uniform grid {2}
Nx cubic cells on x direction, pseudo-periodic condition, non-uniform grid {2}
Nz cubic cells on z direction, opened on top and bottom {2}

Ex dx dt surfaces, -Dx dz lines
1:(Nt+1)*(Nx+1)*Nz

By dz dx surfaces, Hy dt lines
(Nt+1)*(Nx+1)*Nz + UnitRange(1,(Nt+1)*(Nx+1)*Nz)

Ez dz dt surfaces, Dz dx lines
2*(Nt+1)*(Nx+1)*Nz + UnitRange(1,(Nt+1)*(Nx+1)*(Nz+1))

Relations constitutives, cas approx 0

 Ex Dx : Nt×Nx×(Nz+1) [-Ntxz - Ntx] {-12}
 By Hy : Nt×Nx×Nz [-Ntxz] {-8}
 Ez Dz : Nt×Nx×Nz [-Ntxz] {-8}

Équations de flux :

 Phim : Nt×Nx×Nz (cubes à 6 termes) [-Ntxz] {-8}
 
 Phie_t : Nt×Nx×(Nz-1) [-Ntxz + Ntx] {-8 + 4 = -4}
 Phie_x : Nt×Nx×(Nz-1) [-Ntxz + Ntx] {-8 + 4 = -4}
 Phie_z : Nt×Nx×Nz [-Ntxz] {-8}
 
Paramètres libres : [3 Ntx] double couche en haut, simple couche en bas




=#

# phases conditions

NT = 1; NTT = 100
NX = 2
NZ = 3

P = 1 # constitutive approximation level

PT = 1.0
PX = 1.0
PZ = 1.5

DT = PT/NTT
DX = PX/NX
DZ = PZ/NZ

KT = 2*pi
KX = 2*pi*0.1

phitp = exp(-1im*KT*DT*NT) # - for time
phitm = 1.0/phitp
phixp = exp(1im*KX*PX)
phixm = 1.0/phixp

assert(abs(phitp^NTT - 1.0) < 1e-14)

NCTE = NT+2*P+2
NCXE = NX+2*P+2
NCZE = NZ+2*P+2
assert(NCTE < 2^15)
assert(NCXE < 2^15)
assert(NCZE < 2^15)

# position of rectangular cell centers 
CCT = Int16
CCX = Int16
CCZ = Int16
t1 = DT*((1:NCTE) - P - 1.5)
x1 = DX*((1:NCXE) - P - 1.5)
z1 = DZ*((1:NCZE) - P - 1.5)
isgoodCCT(ii::Int) = 1 <= ii <= NCTE
isgoodCCX(ii::Int) = 1 <= ii <= NCXE
isgoodCCZ(ii::Int) = 1 <= ii <= NCZE

# position of segment centers (cubic center faces) 
SCT = Int16
SCX = Int16
SCZ = Int16
t2 = DT*((1:NCTE-1) - P - 1.0)
x2 = DX*((1:NCXE-1) - P - 1.0)
z2 = DZ*((1:NCZE-1) - P - 1.0)
isgoodSCT(ii::Int) = 1 <= ii <= NCTE-1
isgoodSCX(ii::Int) = 1 <= ii <= NCXE-1
isgoodSCZ(ii::Int) = 1 <= ii <= NCZE-1

#=
type ICoord
    it::Int16
    ix::Int16
    iz::Int16
end
=#


# Tuples useful for ind2sub and sub2ind computations
# cell array indices
indsc = (NCTE,NCXE,NCZE)
# segment arrays indices
indst = (NCTE-1,NCXE-1,NCZE-2)
indsx = (NCTE-1,NCXE-1,NCZE-2)
indsz = (NCTE-1,NCXE-1,NCZE-1)

# example:  z1[ind2sub(indsc,sub2ind(indsc,4,4,4))[3]] == 1.25

MAXivte =           (NCZE-2)*(NCXE-1)*(NCTE-1)
MAXivxe = MAXivte + (NCZE-2)*(NCXE-1)*(NCTE-1)
MAXivze = MAXivxe + (NCZE-1)*(NCXE-1)*(NCTE-1)
"""
Return the nature and the index of the segment
"""
function ind2natind(ii::Int)
    if ii <= 0
        throw(BoundsError("indnat", ii))
    elseif ii <= MAXivte
        return (:VTE, ii)
    elseif ii <= MAXivxe
        return (:VXE, ii - MAXivte)
    elseif ii <= MAXivze
        return (:VZE, ii - MAXivxe)
    else
        throw(BoundsError("indnat", ii))
    end
end

function natind2ind(nat, ii)
    if nat == :VTE
        return ii
    elseif nat == :VXE
        return ii + MAXivte
    elseif nat == :VZE
        return ii + MAXivxe
    else
        throw(ErrorException("Bad nat :" * string(nat)))
    end
end

# main vector
vind = Vector{Int64}(NT*NX*(NZ+2*P) + NT*NX*(NZ+2*P) + NT*NX*(NZ+2*P+1))
vindi = zeros(Int64, MAXivze)

ii = 0
for iz in 1:NZ+2*P
    for ix in 1:NX
        for it in 1:NT
            ii += 1
            i0 = natind2ind(:VTE, sub2ind(indst, it+P, ix+P, iz))
            vind[ii] = i0
            vindi[i0] = ii
        end
    end
end
for iz in 1:NZ+2*P
    for ix in 1:NX
        for it in 1:NT
            ii += 1
            i0 = natind2ind(:VXE, sub2ind(indsx, it+P, ix+P, iz))
            vind[ii] = i0
            vindi[i0] = ii
        end
    end
end
for iz in 1:NZ+2*P+1
    for ix in 1:NX
        for it in 1:NT
            ii += 1
            println((it, ix, iz))
            i0 = natind2ind(:VZE, sub2ind(indsz, it+P, ix+P, iz))
            vind[ii] = i0
            vindi[i0] = ii
        end
    end
end

vind0 = Vector{Int64}(MAXivze)
vphi = Vector{Complex{Float64}}(MAXivze)
# Matrices de phase
ii = 0
for iz in 1:NCZE-2
    for ix in 1:NCXE-1
        ix_ = ix - 1 - P
        px, ix__ = fldmod(ix_, NX)
        ix0 = ix__ + P + 1
        for it in 1:NCTE-1
            ii += 1
            it_ = it - 1 - P
            pt, it__ = fldmod(it_, NT)
            it0 = it__ + P + 1
            ind = natind2ind(:VTE, sub2ind(indst, it, ix, iz))
            ind0 = natind2ind(:VTE, sub2ind(indst, it0, ix0, iz))
            println(((it, it0), (ix, ix0), iz, ind0, vindi[ind0]))
            vind0[ii] = vindi[ind0]
            vphi[ii] = phitp^pt * phixp^px
        end
    end
end
for iz in 1:NCZE-2
    for ix in 1:NCXE-1
        ix_ = ix - 1 - P
        px, ix__ = fldmod(ix_, NX)
        ix0 = ix__ + P + 1
        for it in 1:NCTE-1
            ii += 1
            it_ = it - 1 - P
            pt, it__ = fldmod(it_, NT)
            it0 = it__ + P + 1
            ind = natind2ind(:VXE, sub2ind(indsx, it, ix, iz))
            ind0 = natind2ind(:VXE, sub2ind(indsx, it0, ix0, iz))
            println((ind0, vindi[ind0]))
            vind0[ii] = vindi[ind0]
            vphi[ii] = phitp^pt * phixp^px
        end
    end
end
for iz in 1:NCZE-1
    for ix in 1:NCXE-1
        ix_ = ix - 1 - P
        px, ix__ = fldmod(ix_, NX)
        ix0 = ix__ + P + 1
        for it in 1:NCTE-1
            ii += 1
            it_ = it - 1 - P
            pt, it__ = fldmod(it_, NT)
            it0 = it__ + P + 1
            ind = natind2ind(:VZE, sub2ind(indsz, it, ix, iz))
            ind0 = natind2ind(:VZE, sub2ind(indsz, it0, ix0, iz))
            println((ind0, vindi[ind0]))
            vind0[ii] = vindi[ind0]
            vphi[ii] = phitp^pt * phixp^px
        end
    end
end

matphi = sparse(1:MAXivze, vind0, vphi)


# matphi checking
if true
    kz = 0.23
    ex = 0.1
    hy = 0.18
    ez = 0.33
    vvv = Vector{Complex{Float64}}(MAXivze)
    vvv0 = Vector{Complex{Float64}}(NT*NX*(NZ+2*P) + NT*NX*(NZ+2*P) + NT*NX*(NZ+2*P+1))
    
    vvv[:] = NaN

    for iz in 1:NCZE-2
        zz = z1[iz+1]
        for ix in 1:NCXE-1
            xx = x1[ix+1]
            for it in 1:NCTE-1
                tt = t2[it]
                println((it,ix,iz, sub2ind(indst,it,ix,iz)))
                vvv[natind2ind(:VTE, sub2ind(indst,it,ix,iz))] = hy*exp(1im*(-KT*tt + KX*xx + kz*zz))
            end
        end
    end
    for iz in 1:NCZE-2
        zz = z1[iz+1]
        for ix in 1:NCXE-1
            xx = x2[ix]
            for it in 1:NCTE-1
                tt = t1[it+1]
                println((it,ix,iz, sub2ind(indsx,it,ix,iz)))
                vvv[natind2ind(:VXE, sub2ind(indsx,it,ix,iz))] = ez*exp(1im*(-KT*tt + KX*xx + kz*zz))
            end
        end
    end
    for iz in 1:NCZE-1
        zz = z2[iz]
        for ix in 1:NCXE-1
            xx = x1[ix+1]
            for it in 1:NCTE-1
                tt = t1[it+1]
                println((it,ix,iz, sub2ind(indsz,it,ix,iz)))
                vvv[natind2ind(:VZE, sub2ind(indsz,it,ix,iz))] = ex*exp(1im*(-KT*tt + KX*xx + kz*zz))
            end
        end
    end

    for ii in 1:length(vvv0)
        vvv0[ii] = vvv[vind[ii]]
    end

    assert(maximum(abs.(vvv - matphi*vvv0)) < 1e-13)
end

