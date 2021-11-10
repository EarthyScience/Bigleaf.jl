@testset "Reynolds_Number" begin
    Tair,pressure,ustar,z0m = 25,100,0.5,0.5
    R = @inferred Reynolds_Number(Tair,pressure,ustar,z0m)                             
    @test ≈(R, 15870, rtol=1e-3) 
end

@testset "roughness_parameters" begin
    zh = thal.zh
    zr = thal.zr
    LAI = thal.LAI
    keys_exp = (:d, :z0m, :z0m_se)
    rp = @inferred roughness_parameters(Val(:canopy_height), zh)
    #round.(values(rp); sigdigits = 4)
    @test keys(rp) == keys_exp
    @test all(isapproxm.(values(rp), (18.55, 2.65, missing), rtol=1e-3))
    #
    rp = @inferred roughness_parameters(Val(:canopy_height_LAI), zh, LAI)
    #round.(values(rp); sigdigits = 4)
    @test keys(rp) == keys_exp
    @test all(isapproxm.(values(rp), (21.77, 1.419, missing), rtol=1e-3))
    #
    df = copy(tha48)
    dfd = disallowmissing(df)
    #df.wind[1] = missing
    d=0.7*zh
    psi_m = (@inferred stability_correction(columntable(dfd), zr, d)).psi_m
    psi_m = stability_correction(df, zr, d).psi_m
    # note: must use columntable for type stability - but needs compilation timede
    rp = @inferred roughness_parameters(Val(:wind_profile), columntable(dfd), zh, zr; psi_m)
    #round.(values(rp); sigdigits = 4)
    @test keys(rp) == keys_exp
    #@test all(isapproxm.(values(rp), (18.55, 1.879, 0.3561), rtol=1e-3))
    #from R:
    @test all(isapproxm.(values(rp), (18.55, 1.879402, 0.356108), rtol=1e-3))
    #
    # no stability correction
    rp0 = @inferred roughness_parameters(Val(:wind_profile), columntable(dfd), zh, zr; 
        stab_formulation = Val(:no_stability_correction))
    @test keys(rp0) == keys_exp
    # same magnitude as with stability correction
    @test all(isapproxm.(values(rp0), values(rp), rtol=0.5))
    #
    # estimate psi
    #@code_warntype stability_correction(columntable(df), zr, 0.7*zh)
    #@code_warntype roughness_parameters(Val(:wind_profile), columntable(df), zh, zr)
    rp_psiauto = @inferred roughness_parameters(Val(:wind_profile), columntable(dfd), zh, zr)
    @test propertynames(df) == propertynames(tha48) # not changed
    @test rp_psiauto == rp
end

@testset "wind_profile" begin
    datetime, ustar, Tair, pressure, H = values(tha48[1,:])
    z = 30
    d=0.7*thal.zh
    z0m= 2.65
    u30 = @inferred wind_profile(z, ustar, d, z0m)
    @test ≈(u30, 1.93, rtol = 1/100 ) # from R
    #
    u30c = @inferred wind_profile(Val(:Dyer_1970), z, ustar, Tair,pressure, H, d, z0m)
    @test ≈(u30c, 2.31, rtol = 1/100 ) # from R
    #
    z0m=1.9 #2.14 #2.65
    u30 = @inferred wind_profile(z, ustar, d, z0m) # used below
    u30c = @inferred wind_profile(Val(:Dyer_1970), z, ustar, Tair,pressure, H, d, z0m)
    df = copy(tha48)
    dfd = disallowmissing(df)
    windz = @inferred wind_profile(columntable(dfd), z, d, z0m; stab_formulation = Val(:no_stability_correction))    
    windz = wind_profile(df, z, d, z0m; stab_formulation = Val(:no_stability_correction))    
    @test length(windz) == 48
    @test windz[1] == u30
    windzc = @inferred wind_profile(columntable(dfd), z, d, z0m; stab_formulation = Val(:Dyer_1970))    
    @test windzc[1] == u30c
    #plot(windz)
    #plot!(windz2)
    psi_m = stability_correction(df, z, d).psi_m
    windzc2 = @inferred wind_profile(columntable(dfd), z, d, z0m, psi_m)    
    @test windzc2 == windzc
    #
    # estimate z0m
    # need to give zh and zr in addition to many variables in df
    @test_throws Exception wind_profile(df, z, d)    
    windzc3 = @inferred wind_profile(columntable(dfd), z, d; zh=thal.zh, zr=thal.zr, 
        stab_formulation = Val(:Dyer_1970))    
    # may have used slightly different estimated z0m
    #windzc3 - windzc
    @test all(isapprox.(windzc3, windzc, atol=0.1))
    @test windzc3[1] ≈ 2.764203 rtol=1e-3 # from R
end

