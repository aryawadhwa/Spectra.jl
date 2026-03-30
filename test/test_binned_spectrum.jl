using Test
using Spectra

# ===========================================================================
# Tests for BinnedSpectrum
# ===========================================================================

@testset "BinnedSpectrum construction" begin
    edges  = [0.3 0.5; 0.5 1.0; 1.0 2.0; 2.0 5.0; 5.0 10.0]
    counts = [120.0, 340.0, 280.0, 95.0, 30.0]

    spec = BinnedSpectrum(edges, counts)
    @test length(spec) == 5
    @test size(spec.spectral_axis) == (5, 2)
end

@testset "BinnedSpectrum bin_centers and bin_widths" begin
    edges  = [1.0 2.0; 2.0 4.0; 4.0 8.0]
    counts = [10.0, 20.0, 30.0]
    spec   = BinnedSpectrum(edges, counts)

    centers = bin_centers(spec)
    @test centers ≈ [1.5, 3.0, 6.0]

    widths = bin_widths(spec)
    @test widths ≈ [1.0, 2.0, 4.0]
end

@testset "BinnedSpectrum count_rate" begin
    edges    = [0.5 1.0; 1.0 2.0]
    counts   = [500.0, 1000.0]
    exposure = 1000.0
    spec     = BinnedSpectrum(edges, counts; exposure=exposure)

    rates = count_rate(spec)
    @test rates ≈ [0.5, 1.0]
end

@testset "BinnedSpectrum invalid construction" begin
    # Bad edges shape
    @test_throws ArgumentError BinnedSpectrum([1.0 2.0 3.0], [1.0], Dict{Symbol,Any}())
    # Inverted bin
    @test_throws ArgumentError BinnedSpectrum([2.0 1.0], [1.0], Dict{Symbol,Any}())
    # Mismatched flux length
    @test_throws ArgumentError BinnedSpectrum([0.5 1.0; 1.0 2.0], [1.0], Dict{Symbol,Any}())
end

# ===========================================================================
# Tests for an in-memory SingleSpectrum (no file I/O needed)
# ===========================================================================

@testset "SingleSpectrum basic arithmetic" begin
    wave = collect(4000.0:10.0:7000.0)
    flux = ones(length(wave)) .* 100.0
    spec = spectrum(wave, flux)

    scaled = spec * 2.0
    @test maximum(flux_axis(scaled)) ≈ 200.0

    shifted = spec + 50.0
    @test minimum(flux_axis(shifted)) ≈ 150.0
end

@testset "SingleSpectrum metadata passthrough" begin
    wave = collect(1.0:1.0:10.0)
    flux = collect(1.0:1.0:10.0)
    spec = spectrum(wave, flux; instrument="ACIS-S", exposure=5000.0)

    @test spec.instrument == "ACIS-S"
    @test spec.exposure   == 5000.0
end
