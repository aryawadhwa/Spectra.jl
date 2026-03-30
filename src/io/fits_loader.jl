# src/io/fits_loader.jl
#
# Prototype FITS data loading routines for Spectra.jl.
# This demonstrates the proposed `Spectra.load()` API and provides
# a concrete starting point for the OGIP migration work.
#
# Dependencies: FITSIO.jl (to be added to Project.toml)
#
# Author: Arya Wadhwa (GSoC 2026 exploration)

using FITSIO

"""
    load(filepath::String; format=:auto) -> AbstractSpectrum

Load an astrophysical spectrum from a file. Automatically detects the format
from the file extension and FITS header keywords, or accepts an explicit
`format` keyword.

Supported formats:
- `:fits1d` — Standard 1D optical/UV spectrum (FITS binary table or image)
- `:ogip`   — X-ray OGIP PHA spectrum (requires associated ARF/RMF files)
- `:auto`   — Attempts automatic detection (default)

# Examples

```julia
spec = Spectra.load("spectrum.fits")
spec = Spectra.load("src.pha", format=:ogip)
```
"""
function load(filepath::String; format::Symbol=:auto)
    if !isfile(filepath)
        error("File not found: $filepath")
    end
    fmt = format == :auto ? _detect_format(filepath) : format
    if fmt == :fits1d
        return _load_fits1d(filepath)
    elseif fmt == :ogip
        return _load_ogip_pha(filepath)
    else
        error("Unsupported format: $fmt. Use :fits1d, :ogip, or :auto.")
    end
end

"""
Internal: Detect spectral format from file extension and FITS header keywords.
"""
function _detect_format(filepath::String)
    ext = lowercase(splitext(filepath)[2])
    ext in (".pha", ".pi") && return :ogip
    # For .fits files, check the HDUCLAS1 keyword
    if ext in (".fits", ".fit", ".fts")
        f = FITS(filepath, "r")
        try
            if length(f) >= 2
                hdr = read_header(f[2])
                if haskey(hdr, "HDUCLAS1") && hdr["HDUCLAS1"] == "SPECTRUM"
                    return :ogip
                end
            end
        finally
            close(f)
        end
        return :fits1d
    end
    return :fits1d
end

"""
Internal: Load a standard 1D optical/UV FITS spectrum.
Expects a binary table HDU with WAVELENGTH (or WAVE) and FLUX columns.
"""
function _load_fits1d(filepath::String)
    f = FITS(filepath, "r")
    meta = Dict{Symbol, Any}()

    # Extract primary header metadata
    primary_hdr = read_header(f[1])
    for key in ("TELESCOP", "INSTRUME", "OBJECT", "DATE-OBS", "EXPTIME")
        if haskey(primary_hdr, key)
            meta[Symbol(lowercase(key))] = primary_hdr[key]
        end
    end

    # Try reading binary table from HDU 2
    if length(f) < 2
        close(f)
        error("No binary table extension found in $filepath. Cannot parse as 1D spectrum.")
    end

    table = f[2]
    cols  = FITSIO.colnames(table)
    close(f)

    wave_col = findfirst(c -> uppercase(c) in ("WAVELENGTH", "WAVE", "LAMBDA"), cols)
    flux_col = findfirst(c -> uppercase(c) in ("FLUX", "FLUX_DENSITY"), cols)

    if isnothing(wave_col) || isnothing(flux_col)
        error("Could not locate WAVELENGTH/FLUX columns. Found: $(join(cols, \", \"))")
    end

    f2 = FITS(filepath, "r")
    wave = read(f2[2], cols[wave_col])
    flux = read(f2[2], cols[flux_col])
    close(f2)

    return spectrum(wave, flux; meta...)
end

"""
Internal: Load an OGIP PHA spectral file.
Reads CHANNEL and COUNTS (or RATE) columns from the SPECTRUM extension.
Does not yet load ARF or RMF — that is Phase 1 of the GSoC project.
"""
function _load_ogip_pha(filepath::String)
    f = FITS(filepath, "r")
    meta = Dict{Symbol, Any}()

    # Find the SPECTRUM extension
    spectrum_hdu = nothing
    for i in 1:length(f)
        hdr = read_header(f[i])
        if haskey(hdr, "HDUCLAS1") && hdr["HDUCLAS1"] == "SPECTRUM"
            spectrum_hdu = i
            # Extract OGIP-standard metadata
            for key in ("TELESCOP", "INSTRUME", "OBJECT", "EXPOSURE", "RESPFILE", "ANCRFILE")
                haskey(hdr, key) && (meta[Symbol(lowercase(key))] = hdr[key])
            end
            break
        end
    end

    if isnothing(spectrum_hdu)
        close(f)
        error("No SPECTRUM extension found in $filepath. Is this a valid OGIP PHA file?")
    end

    cols    = FITSIO.colnames(f[spectrum_hdu])
    channel = read(f[spectrum_hdu], "CHANNEL")

    # PHA files store either COUNTS (integer) or RATE (float)
    flux_col = "COUNTS" in cols ? "COUNTS" : "RATE"
    flux     = Float64.(read(f[spectrum_hdu], flux_col))
    close(f)

    meta[:flux_type] = flux_col

    # NOTE: For a proper OGIP loader, channels must be mapped to energies
    # using the RMF (Redistribution Matrix File). This is a stub that
    # returns the raw channel-space spectrum. Full energy calibration
    # is the primary Phase 1 goal of this GSoC project.
    return spectrum(Float64.(channel), flux; meta...)
end
