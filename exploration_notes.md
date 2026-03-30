# Spectra.jl: Codebase Exploration Notes

**Author:** Arya Wadhwa  
**Date:** March 2026  
**Context:** GSoC 2026 application exploration for the "Spectra.jl across the electromagnetic spectrum" project.

---

## 1. Repository Structure

```
src/
  Spectra.jl              # Module root: abstract types, arithmetic, dispatch
  spectrum_single.jl      # SingleSpectrum concrete type
  spectrum_echelle.jl     # EchelleSpectrum concrete type (multi-order)
  spectrum_ifu.jl         # IFUSpectrum concrete type (3D integral field)
  utils.jl                # Blackbody and utility functions
  plotting.jl             # RecipesBase plotting integration
  transforms/
    resampler.jl          # SpectrumResampler
    transforms.jl         # redden / deredden (dust extinction)
  old_functionality/
    fitting/              # Deprecated fitting routines (commented out in module)
```

---

## 2. Core Data Structures

### AbstractSpectrum{S, F}

The root abstract type is parameterized over two numeric types:
- `S` — the element type of the spectral axis (e.g. wavelength in Angstroms, or energy in keV)
- `F` — the element type of the flux axis

Every concrete subtype **must** implement three fields:
```julia
spectral_axis::Array{S, M}
flux_axis::Array{F, N}
meta::Dict{Symbol, Any}
```

The `meta` dictionary is clever: it is transparently accessible as properties via `Base.getproperty` overloading, so users can write `spec.name` instead of `spec.meta[:name]`. This is a clean design choice.

### Spectrum{S, F, M, N}

The primary concrete type is parameterized over four type parameters: element types S, F and array dimensions M, N. A single constructor dispatches to `SingleSpectrum`, `EchelleSpectrum`, or `IFUSpectrum` based on the shape of the inputs.

**Validation on construction:**
- Dimensional size check: `size(s, 1) == size(f, 1)` is enforced.
- Monotonicity check: The spectral axis must be strictly increasing or decreasing across all columns (checked via `eachcol`).

**Observation:** The TODO comment inside the constructor mentions investigating Holy Traits for validation. This is a meaningful future direction as the number of spectrum types grows.

### Existing Spectrum Types

| Type | spectral_axis | flux_axis | Use Case |
|---|---|---|---|
| `SingleSpectrum` | M-vector | M-vector | Standard 1D optical/UV/IR/X-ray |
| `EchelleSpectrum` | M x N matrix | M x N matrix | Multi-order echelle gratings |
| `IFUSpectrum` | M-vector | M x N x K cube | Integral field spectrographs |

A fourth type, `BinnedSpectrum`, is referenced in comments but **not yet implemented**. This is the most relevant type for X-ray OGIP spectra, where the spectral axis is defined by energy bin edges rather than center points.

---

## 3. Unit Handling

The package integrates `Unitful.jl` throughout. The `spectrum()` constructor asserts that a `Quantity`-typed spectral axis must carry dimensions of either wavelength (𝐋) or energy (𝐋² 𝐌 𝐓⁻²). This covers both UV/optical (Angstroms) and X-ray (keV) regimes correctly.

`UnitfulAstro.jl` is used in documentation examples. However, the current unit conversion engine appears minimal — `ustrip` removes units, but there is no in-place conversion between common astronomical flux notations (e.g., Jansky to erg/s/cm²/Hz). This is one of the things the GSoC project should add.

---

## 4. The Missing OGIP Data Loader

The most significant gap visible from this exploration: **there is no data loader anywhere in the package.** There is no `load()`, `read_fits()`, or `from_ogip()` function. The package can model and manipulate spectra that are already in memory, but it has no ingestion pipeline.

In comparison, `SpectralFitting.jl` already contains a working OGIP parser. The migration target is roughly:

```
SpectralFitting.jl/src/datasets/ogip/
  -> Spectra.jl/src/io/ogip.jl  (proposed)
```

The OGIP standard (OGIP/92-007) defines PHA files storing counts per channel, ARF files storing effective area, and RMF files storing the energy redistribution matrix. A correct X-ray spectral loader must handle all three simultaneously.

---

## 5. Open Questions for Mentors

1. **BinnedSpectrum:** Is `BinnedSpectrum` the correct target type for X-ray OGIP ingested data, given that OGIP spectra have bin edges, not center points? Or should the existing `SingleSpectrum` be extended?
2. **IO module location:** Should all loaders live under a new `src/io/` directory with one file per format, or should they be integrated directly into the respective spectrum type files?
3. **Dependency budget:** Is adding `FITSIO.jl` as a direct dependency acceptable? Or should FITS loading be placed in a package extension to keep the core lightweight?
4. **Test data library:** The issue mentions a curated library of spectral data files. Where is this hosted so TDD suites can be written against real instrument data?

---

## 6. Proposed First Contribution

Based on this exploration, the most actionable starting point is to propose the `src/io/` directory structure and implement a minimal FITS reader for a standard 1D optical spectrum (the simplest format). This proves the ingestion pattern before tackling the complexity of OGIP response matrices.

A minimal interface could look like:

```julia
# Proposed public API
spec = Spectra.load("spectrum.fits")          # auto-detect format
spec = Spectra.load("src.pha", format=:ogip)  # explicit format
```
