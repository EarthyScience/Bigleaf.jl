abstract type EsatMethod end
struct Sonntag1990 <: EsatMethod end
struct Alduchov1996 <: EsatMethod end
struct Allen1998 <: EsatMethod end

"""
    Esat_slope(Tair; Esat_formula, constants) 
    Esat_from_Tair(Tair; Esat_formula, constants) 
    Esat_from_Tair_deriv(Tair; Esat_formula, constants) 

Saturation Vapor Pressure (Esat) and Slope of the Esat Curve

# Arguments
- `Tair`:      Air temperature (deg C)
- `Esat_formula=Sonntag1990()`:   Esat_formula to be used. Either 
   `Sonntag1990()` (Default), `Alduchov1996()`, or `Allen1998()`.
- `constants=BigleafConstants()`: Dictionary with entry  :Pa2kPa 

# Details
Esat (kPa) is calculated using the Magnus equation:

```Esat = a * exp((b * Tair) / (c + Tair)) / 1000```
where the coefficients a, b, c take different values depending on the Esat_formula use
The default values are from Sonntag 1990 (a=611.2, b=17.62, c=243.12). This versi
of the Magnus equation is recommended by the WMO (WMO 2008; p1.4-29). Alternativel
parameter values determined by Alduchov & Eskridge 1996 or Allen et al. 1998 can b
used (see references).
The slope of the Esat curve (``\\Delta``) is calculated as the first derivative of the function:
 
``\\Delta = {dEsat \\over dTair}``

# Value   
- `Esat_from_Tair`: Saturation vapor pressure (kPa)
- `Esat_from_Tair_deriv`: Slope of the saturation vapor pressure curve (kPa K-1)
- `Esat_slope`: A tuple of the two above values
# References
Sonntag D. 1990: Important new values of the physical constants of 1986, vapor 
pressure formulations based on the ITS-90 and psychrometric formulae. 
Zeitschrift fuer Meteorologie 70, 340-344.

World Meteorological Organization 2008: Guide to Meteorological Instruments
and Methods of Observation (WMO-No.8). World Meteorological Organization,
Geneva. 7th Edition,

Alduchov, O. A. & Eskridge, R. E., 1996: Improved Magnus form approximation of 
saturation vapor pressure. Journal of Applied Meteorology, 35, 601-609

Allen, R.G., Pereira, L.S., Raes, D., Smith, M., 1998: Crop evapotranspiration -
Guidelines for computing crop water requirements - FAO irrigation and drainage
paper 56, FAO, Rome.

```@example
Esat_from_Tair(20.0)          # Esat in kPa
Esat_from_Tair_deriv(20.0)    # its derivative to temperature in kPa K-1
Esat_slope(20.0)              # both as a tuple
```
"""
function Esat_slope(Tair::Number; Esat_formula=Sonntag1990(), constants=BigleafConstants()) 
  Esat = Esat_from_Tair(Tair; Esat_formula, constants)
  Delta = Esat_from_Tair_deriv(Tair; Esat_formula, constants)
  Esat, Delta
end,
function Esat_from_Tair(Tair::FT; Esat_formula=Sonntag1990(), constants=BigleafConstants()) where FT
  a,b,c = map(x -> convert(FT,x), get_EsatCoef(Esat_formula))
  Esat = a * exp((b * Tair) / (c + Tair)) * FT(constants.Pa2kPa)
end,
function Esat_from_Tair_deriv(Tair::FT; Esat_formula=Sonntag1990(), constants=BigleafConstants()) where FT
  # slope of the saturation vapor pressure curve
  #Delta = eval(D(expression(a * exp((b * Tair) / (c + Tair))),name="Tair"))
  a,b,c = map(x -> convert(FT,x), get_EsatCoef(Esat_formula))
  #Delta_Pa = @. a*(b / (Tair + c) + (-Tair*b) / ((Tair + c)^2))*exp((Tair*b) / (Tair + c))
  Delta_Pa = a * (exp((b * Tair)/(c + Tair)) * (b/(c + Tair) - (b * Tair)/(c + Tair)^2))
  Delta = Delta_Pa .* convert(FT, constants.Pa2kPa)
end
    
get_EsatCoef(::Sonntag1990) = (a=611.2,b=17.62,c=243.12)
get_EsatCoef(::Alduchov1996) = (a=610.94,b=17.625,c=243.04)
get_EsatCoef(::Allen1998) = (a=610.8,b=17.27,c=237.3)

"""
    LE_to_ET(LE,Tair)
    ET_to_LE(ET,Tair)

Convert evaporative water flux from mass (ET=evapotranspiration)
             to energy (LE=latent heat flux) units, or vice versa.

# Arguments
- LE   Latent heat flux (W m-2)
- ET   Evapotranspiration (kg m-2 s-1)
- Tair Air temperature (deg C)

# Details
The conversions are given by:
- ``ET = LE/\\lambda``
- ``LE = \\lambda ET``

where ``\\lambda`` is the latent heat of vaporization (J kg-1) as calculated by
[`latent_heat_vaporization`](@ref).

# Examples
```jldoctest; output = false
# LE of 200 Wm-2 and air temperature of 25degC
ET = LE_to_ET(200,25)
≈(ET, 8.19e-5, atol =1e-7)
# output
true
```
"""
function LE_to_ET(LE,Tair)
  lambda = latent_heat_vaporization(Tair)
  ET     = LE/lambda
end,
function ET_to_LE(ET,Tair)
  lambda = latent_heat_vaporization(Tair)
  LE     = ET*lambda
end


"""
    ms_to_mol(G_ms,Tair,pressure; constants=BigleafConstants())
    mol_to_ms(G_mol,Tair,pressure; constants=BigleafConstants())

Converts conductances from mass (m s-1) to molar units (mol m-2 s-1), or vice versa

# Details
The conversions are given by
- ``G_{mol} = G_{ms}  \\, pressure / (Rgas  Tair)``
- ``G_{ms} = G_{mol}  \\, (Rgas  Tair) / pressure``

where Tair is in deg C and pressure in Pa (converted from kPa internally).

# References 
Jones, H_G_ 1992_ Plants and microclimate: a quantitative approach to 
  environmental plant physiology_
  2nd Edition_, Cambridge University Press, Cambridge_ 428 

# Examples
```jldoctest; output = false
G_ms,Tair,pressure = 0.005,25.0,100.0
rmol = ms_to_mol(G_ms,Tair,pressure)
≈(rmol, 0.2017, atol =1e-4)
# output
true
```
"""
function ms_to_mol(G_ms,Tair,pressure; constants=BigleafConstants())
  Tair     = Tair + oftype(Tair,constants.Kelvin)
  pressure = pressure * oftype(Tair, constants.kPa2Pa)
  G_mol  = G_ms * pressure / (oftype(Tair, constants.Rgas)  * Tair)
end,
function mol_to_ms(G_mol,Tair,pressure; constants=BigleafConstants())
  Tair     = Tair + oftype(Tair,constants.Kelvin)
  pressure = pressure * oftype(Tair, constants.kPa2Pa)
  G_ms  = G_mol * (oftype(Tair, constants.Rgas) * Tair) / (pressure)
end


"""
    VPD_to_rH(VPD,Tair; ...)
    H_to_VPD(rH,Tair; ...)
    e_to_rH(e,Tair; ...)
    VPD_to_e(VPD,Tair; ...)
    e_to_VPD(e,Tair; ...)
    e_to_q(e,pressure; ...)
    q_to_e(q,pressure; ...)
    q_to_VPD(q,Tair,pressure; ...)
    VPD_to_q(VPD,Tair,pressure; ...)

Conversion between vapor pressure (e), vapor pressure deficit (VPD),
             specific humidity (q), and relative humidity (rH).

# Arguments             
- Tair:      Air temperature (deg C)
- pressure:  Atmospheric pressure (kPa)
- e:         Vapor pressure (kPa)
- q:        Specific humidity (kg kg-1)
- VPD:       Vapor pressure deficit (kPa)
- rH:        Relative humidity (-)

All functions accept the optional arguments:
- `Esat_formula`: formula used in [`Esat_from_Tair`](@ref)
- `constants`: dictionary from [`BigleafConstants`](@ref) with entries
  eps and Pa2kPa

# Rreferences 
Foken, T, 2008: Micrometeorology_ Springer, Berlin, Germany.
"""
function VPD_to_rH(VPD,Tair; Esat_formula=Sonntag1990(),
                      constants=BigleafConstants())
  esat = Esat_from_Tair(Tair; Esat_formula = Esat_formula,constants)
  rH   = 1 - VPD/esat
end,
function rH_to_VPD(rH,Tair; Esat_formula=Sonntag1990(),
                      constants=BigleafConstants())
  if !ismissing(rH) && rH > 1 
    @warn("Expected relative humidity (rH) between 0 and 1, but was ", rH)
  end
  esat = Esat_from_Tair(Tair; Esat_formula = Esat_formula,constants)
  VPD  = esat - rH*esat
end,
function e_to_rH(e,Tair; Esat_formula=Sonntag1990(),
                    constants=BigleafConstants())
  esat = Esat_from_Tair(Tair; Esat_formula = Esat_formula,constants)
  if !ismissing(e) && (e > esat + sqrt(eps()))
    @warn("Provided vapour pressure that was higher than saturation_
             Returning rH=1 for those cases.")
  end
  rH  = min(1, e/esat)
end,
function VPD_to_e(VPD,Tair; Esat_formula=Sonntag1990(),
                     constants=BigleafConstants())
  esat = Esat_from_Tair(Tair; Esat_formula = Esat_formula,constants)
  e    = esat - VPD
end,
function e_to_VPD(e,Tair; Esat_formula=Sonntag1990(),
                     constants=BigleafConstants())
  esat = Esat_from_Tair(Tair; Esat_formula = Esat_formula,constants)
  VPD  = esat - e
end,
function e_to_q(e::FT,pressure; constants=BigleafConstants()) where FT
  q = (FT(constants.eps) * e / (pressure - (FT(1)-FT(constants.eps)) * e))
end,
function q_to_e(q::FT,pressure; constants=BigleafConstants()) where FT
  e = (q * pressure / ((FT(1)-FT(constants.eps)) * q + FT(constants.eps)))
end,
function q_to_VPD(q,Tair,pressure; Esat_formula=Sonntag1990(),
                     constants=BigleafConstants())
  esat = Esat_from_Tair(Tair; Esat_formula = Esat_formula,constants)
  e    = q_to_e(q,pressure; constants)
  VPD  = esat - e
end,
function VPD_to_q(VPD,Tair,pressure; Esat_formula=Sonntag1990(),
                     constants=BigleafConstants())
  esat = Esat_from_Tair(Tair; Esat_formula = Esat_formula,constants)
  e    = esat - VPD
  q    = e_to_q(e,pressure;constants)
end

"""
    Rg_to_PPFD(Rg,J_to_mol=4.6,frac_PAR=0.5)
    PPFD_to_Rg(PPFD,J_to_mol=4.6,frac_PAR=0.5)

Conversions between Global Radiation (W m-2) and Photosynthetic Photon Flux 
Density (umol m-2 s-1)

# Arguments
- Rg:       Global radiation = incoming short-wave radiation at the surface (W m-2)
- PPFD:     Photosynthetic photon flux density (umol m-2 s-1)
- J_to_mol: Conversion factor from J m-2 s-1 (= W m-2) to umol (quanta) m-2 s-1
- frac_PAR: Fraction of incoming solar irradiance that is photosynthetical-          
  active radiation (PAR); defaults to 0.5

# Details
The conversion is given by:

 ``PPFD = Rg * frac_PAR * J_to_mol``

by default, the combined conversion factor (`frac_PAR * J_to_mol`) is 2.3

# Examples
```@example
# convert a measured incoming short-wave radiation of 500 Wm-2 to
# PPFD in umol m-2 s-1 and backwards
Rg_to_PPFD(500)
PPFD_to_Rg(1150)
```
"""
function Rg_to_PPFD(Rg,J_to_mol=4.6,frac_PAR=0.5)
  PPFD = Rg * frac_PAR * J_to_mol
end,
function PPFD_to_Rg(PPFD,J_to_mol=4.6,frac_PAR=0.5)
  Rg = PPFD / frac_PAR / J_to_mol
end

"""
    kg_to_mol(mass, molarMass=BigleafConstants()[:H2Omol])

Conversion between Mass (kg) and Molar Units (mol).
"""
function kg_to_mol(mass, molarMass=BigleafConstants()[:H2Omol])
  moles = mass / molarMass
end

"""
    umolCO2_to_gC(CO2_flux; constants=BigleafConstants())
    gC_to_umolCO2(C_flux; constants=BigleafConstants())


Convert CO2 quantities from (umol CO2 m-2 s-1) to (g C m-2 d-1) and vice versa.

# Arguments
- CO2_flux  CO2 flux (umol CO2 m-2 s-1)
- C_flux    Carbon (C) flux (gC m-2 d-1)
- `constants`: dictionary from [`BigleafConstants`](@ref) with entries:
  Cmol, umol2mol, mol2umol, kg2g, g2kg, says2seconds

# Examples                 
```@example
umolCO2_to_gC(20)  # gC m-2 d-1
```
"""
function umolCO2_to_gC(CO2_flux::FT; constants=BigleafConstants()) where FT
  C_flux = CO2_flux * FT(constants.umol2mol) * FT(constants.Cmol) * 
  FT(constants.kg2g) * FT(constants.days2seconds)
end,
function gC_to_umolCO2(C_flux::FT; constants=BigleafConstants()) where FT
  CO2_flux = (C_flux * FT(constants.g2kg)  / FT(constants.days2seconds)) / 
  FT(constants.Cmol)  * FT(constants.mol2umol)
end
