using Bigleaf, Test

@testset "Esat_from_Tair" begin
  Tair = 15.0
  Esat_formula=Sonntag1990()
  constants=BigleafConstants()
  eSat = Esat_from_Tair(Tair; Esat_formula, constants)
  # regression test
  @test isapprox(eSat, 1.702, atol=1e-3)
  eSat2 = Esat_from_Tair(Tair; Esat_formula = Alduchov1996(), constants)
  eSat3 = Esat_from_Tair(Tair; Esat_formula = Allen1998(), constants)
  # different results with different methods
  @test !(eSat2 ≈ eSat)
  @test !(eSat3 ≈ eSat)
  @test !(eSat3 ≈ eSat2)
  # but very similar
  @test ≈(eSat2, eSat, rtol = 1e-2)
  @test ≈(eSat3, eSat, rtol = 1e-2)
end

@testset "Esat_from_Tair_slope" begin
    step = 0.1
    Tair = range(20.0, 22.0, step = step)
    Esat_formula=Sonntag1990()
    constants=BigleafConstants()
    Esat = Esat_from_Tair.(Tair)
    delta = Esat_from_Tair_deriv.(Tair)
    delta2 = diff(Esat)/step
    @test all(isapprox.(delta[2:end] - delta2, 0, atol=1e-3))
    Esat3, delta3 = Esat_slope(Tair[1])
    @test Esat3 ≈ Esat[1]
    @test delta3 ≈ delta[1]
end

@testset "LE_to_ET" begin
  ET0 = 200.0
  Tair = 25.0
  res = LE_to_ET(ET0, Tair)
  # regression test
  @test ≈(res, 8.19e-5, atol =1e-7)
  ET = ET_to_LE(res, Tair)
  @test ET ≈ ET0
end

@testset "LE_to_ET" begin
  G_ms,Tair,pressure = 0.005,25.0,100.0
  rmol = ms_to_mol(G_ms,Tair,pressure)
  # regression test
  @test ≈(rmol, 0.2017, atol =1e-4)
  rms = mol_to_ms(rmol, Tair, pressure)
  @test rms ≈ G_ms
end



@testset "e_to_rH" begin
  Tair = 25.0
  eSat = Esat_from_Tair(Tair)
  e = range(0,eSat, length = 5)
  vpd = e_to_VPD.(e, Tair)
  @test VPD_to_e.(vpd, Tair) ≈ e
  rH0 = VPD_to_rH.(vpd,Tair)
  @test rH_to_VPD.(rH0, Tair) ≈ vpd
  rH = e_to_rH.(e, Tair)
  @test  rH ≈ rH0
  #
  rHOversat = @test_logs (:warn,r"higher than saturation") e_to_rH(eSat + 1e-3, Tair)
  @test rHOversat == 1.0
  #
  vpd2 = @test_logs (:warn,r"between 0 and 1") rH_to_VPD(1.2, Tair)
  @test vpd2 < 0 # garbage in garbage out
end

@testset "e_to_q" begin
  pressure = 0.1 #kPa
  Tair = 25.0
  eSat = Esat_from_Tair(Tair)
  e = range(0,eSat, length = 5)
  q = e_to_q.(e, pressure)
  @test q_to_e.(q, pressure) ≈ e
  vpd = q_to_VPD.(q, Tair, pressure)
  @test VPD_to_q.(vpd, Tair, pressure) ≈ q
end

@testset "Rg_to_PPFD" begin
  Rg = 500.0
  PPFD = Rg_to_PPFD(Rg)
  @test PPFD_to_Rg(PPFD) ≈ Rg
end

@testset "kg_to_mol" begin
  mass = 10.0
  molarMass =BigleafConstants().H2Omol
  amountOfSubstance = kg_to_mol(mass, molarMass)
  @test amountOfSubstance == mass/molarMass
end

@testset "umolCO2_to_gC" begin
  CO2_flux = 50.0
  C_flux = umolCO2_to_gC(CO2_flux)
  @test gC_to_umolCO2(C_flux) ≈ CO2_flux
end
