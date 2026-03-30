# BinnedSpectrum
#
# This file implements the BinnedSpectrum concrete type, which is referenced
# as a TODO in the main Spectra.jl module. BinnedSpectrum is the correct
# structure for X-ray OGIP spectra, where the spectral axis consists of
# energy bin edges (an M x 2 matrix) rather than central values.
#
# OGIP standard reference: OGIP/92-007 (PHA spectral format)

"""
    BinnedSpectrum <: AbstractSpectrum

A spectrum where the spectral axis defines **bin edges** rather than central
wavelength or energy values. This is the standard representation for X-ray
spectra from observatories such as Chandra, XMM-Newton, and NuSTAR, which
report photon counts accumulated over discrete instrument channels.

The `spectral_axis` is an `M x 2` matrix where each row `[E_lo, E_hi]`
defines the lower and upper boundary of one energy bin. The `flux_axis` is
an `M`-length vector of counts (or count rates) per bin.

# Fields
- `spectral_axis`: `M x 2` matrix of bin edges
- `flux_axis`: `M`-length vector of bin values
- `meta`: metadata dictionary (instrument, exposure, units, etc.)

# Examples

```julia
# Define 5 energy bins from 0.3 to 10 keV
edges = [0.3 0.5; 0.5 1.0; 1.0 2.0; 2.0 5.0; 5.0 10.0]
counts = [120.0, 340.0, 280.0, 95.0, 30.0]
spec = BinnedSpectrum(edges, counts, Dict(:instrument => "ACIS-S"))

# Access bin centers (useful for plotting)
centers = bin_centers(spec)
widths  = bin_widths(spec)
```
"""
mutable struct BinnedSpectrum{S<:Number, F<:Number} <: AbstractSpectrum{S, F}
    spectral_axis::Matrix{S}   # M x 2: each row is [lo, hi] bin edge
    flux_axis::Vector{F}       # M: one value per bin
    meta::Dict{Symbol, Any}

    function BinnedSpectrum{S, F}(edges::Matrix{S}, flux::Vector{F}, meta::Dict{Symbol, Any}) where {S<:Number, F<:Number}
        nbins = size(edges, 1)
        if size(edges, 2) != 2
            throw(ArgumentError("spectral_axis must be an M x 2 matrix of bin edges [lo, hi]."))
        end
        if length(flux) != nbins
            throw(ArgumentError("flux_axis length ($(length(flux))) must equal the number of bins ($nbins)."))
        end
        # Validate all bins are positively oriented
        if any(edges[:, 1] .>= edges[:, 2])
            throw(ArgumentError("All bin edges must satisfy lo < hi."))
        end
        return new{S, F}(edges, flux, meta)
    end
end

function BinnedSpectrum(edges::Matrix{S}, flux::Vector{F}, meta::Dict{Symbol, Any}) where {S<:Number, F<:Number}
    BinnedSpectrum{S, F}(edges, flux, meta)
end

function BinnedSpectrum(edges::Matrix{S}, flux::Vector{F}; kwds...) where {S<:Number, F<:Number}
    BinnedSpectrum{S, F}(edges, flux, Dict{Symbol, Any}(kwds))
end

# Human readable display
function Base.show(io::IO, spec::BinnedSpectrum)
    nbins = size(spec.spectral_axis, 1)
    lo = spec.spectral_axis[1, 1]
    hi = spec.spectral_axis[end, 2]
    println(io, "BinnedSpectrum($(eltype(spec.spectral_axis)), $(eltype(spec.flux_axis)))")
    println(io, "  bins: $nbins  |  range: $lo .. $hi")
    println(io, "  flux axis: $(minimum(spec.flux_axis)) .. $(maximum(spec.flux_axis))")
    print(io, "  meta: ", spec.meta)
end

# Standard array interface
Base.length(spec::BinnedSpectrum)   = length(spec.flux_axis)
Base.size(spec::BinnedSpectrum)     = size(spec.flux_axis)

"""
    bin_centers(spec::BinnedSpectrum)

Return the midpoint of each energy or wavelength bin.
"""
bin_centers(spec::BinnedSpectrum) = (spec.spectral_axis[:, 1] .+ spec.spectral_axis[:, 2]) ./ 2

"""
    bin_widths(spec::BinnedSpectrum)

Return the width of each energy or wavelength bin.
"""
bin_widths(spec::BinnedSpectrum) = spec.spectral_axis[:, 2] .- spec.spectral_axis[:, 1]

"""
    count_rate(spec::BinnedSpectrum)

If `spec.meta[:exposure]` is set (in seconds), return counts per second per bin.
Otherwise throws an error indicating the metadata is missing.
"""
function count_rate(spec::BinnedSpectrum)
    if !haskey(spec.meta, :exposure)
        error("Exposure time not found in spectrum metadata. Set meta[:exposure] in seconds.")
    end
    return spec.flux_axis ./ spec.meta[:exposure]
end
