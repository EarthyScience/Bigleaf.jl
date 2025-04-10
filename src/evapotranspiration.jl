abstract type ETMethod end
struct PriestleyTaylor <: ETMethod end
struct PenmanMonteith <: ETMethod end

"""
    potential_ET(::PriestleyTaylor, Tair, pressure, Rn; G=0.0, S=0.0, ...)
    potential_ET(::PriestleyTaylor, Tair, pressure, Rn, G, S; ...)

    potential_ET(::PenmanMonteith, Tair, pressure, Rn, VPD, Ga_h;
      G=zero(Tair),S=zero(Tair), ...)
    fpotential_ET(::PenmanMonteith, Tair, pressure, Rn, VPD, Ga_h, G, S; ...)
    
    potential_ET!(df, approach; ...)

Potential evapotranspiration according to Priestley & Taylor 1972 or
the Penman-Monteith equation with a prescribed surface conductance.

# Arguments             
- `Tair`:      Air temperature (degC)
- `pressure`:  Atmospheric pressure (kPa)
- `Rn`:        Net radiation (W m-2)
- `VPD`:       Vapor pressure deficit (kPa)
- `Ga`:        Aerodynamic conductance to heat/water vapor (m s-1)
- `df`:        DataFrame with the above variables
- `approach`:  Approach used: Either `PriestleyTaylor()`  or `PenmanMonteith()`.
optional:
- `G=0.0`:         Ground heat flux (W m-2). Defaults to zero.
- `S=0.0`:         Sum of all storage fluxes (W m-2) . Defaults to zero.
- `constants=`[`BigleafConstants`](@ref)`()`: physical constants (cp, eps, Rd, Rgas)
for `PriestleyTaylor`:
- `alpha = 1.26`:   Priestley-Taylor coefficient
for PenmanMonteith:
- `Gs_pot = 0.6`:    Potential/maximum surface conductance (mol m-2 s-1); 
- `Esat_formula`: formula used in [`Esat_from_Tair`](@ref)
  

# Details
Potential evapotranspiration is calculated according to Priestley & Taylor, 1972
(`approach = PriestleyTaylor()`:

``\\mathit{LE}_{pot} = (\\alpha \\, \\Delta \\, (Rn - G - S)) / (\\Delta + \\gamma)``

``\\alpha`` is the Priestley-Taylor coefficient, ``\\Delta`` is the slope 
of the saturation vapor pressure curve (kPa K-1), and ``\\gamma`` is the 
psychrometric constant (kPa K-1).

If `approach = PenmanMonteith()`, potential evapotranspiration is calculated according
to the Penman-Monteith equation:

``\\mathit{LE}_{pot} = (\\Delta \\, (R_n - G - S) + \\rho \\, c_p \\, \\mathit{VPD} \\; G_a) / 
(\\Delta + \\gamma \\, (1 + G_a/G_{s pot})``

where ``\\Delta`` is the slope of the saturation vapor pressure curve (kPa K-1),
``\\rho`` is the air density (kg m-3), 
and ``\\gamma`` is the psychrometric constant (kPa K-1).
The value of ``G_{s pot}`` is typically a maximum value of ``G_s`` observed at the site, 
e.g. the ``90^{th}`` percentile of ``G_s`` within the growing season.

Ground heat flux and storage heat flux `G` or `S` are provided as optional 
arguments. In the input-explicit variants, they default to zero.
In the data-frame arguments, they default to missing, which results
in assuming them to be zero which is displayed in a log-message.
Note that in difference to the bigleaf R package, you explicitly need to
care for missing values (see examples).
      
# Value
NamedTuple with the following entries:
- `ET_pot`: Potential evapotranspiration (kg m-2 s-1)
- `LE_pot`: Potential latent heat flux (W m-2)

# References 
- Priestley, C_H_B., Taylor, R_J., 1972: On the assessment of surface heat flux
  and evaporation using large-scale parameters. Monthly Weather Review 100, 81-92.  
- Allen, R_G., Pereira L_S., Raes D., Smith M., 1998: Crop evapotranspiration -
  Guidelines for computing crop water requirements - FAO Irrigation and drainage
  paper 56.
- Novick, K_A., et al. 2016: The increasing importance of atmospheric demand
  for ecosystem water and carbon fluxes. Nature Climate Change 6, 1023 - 1027.
            
# See also 
[`surface_conductance`](@ref)

# Examples
Calculate potential ET of a surface that receives a net radiation of 500 Wm-2
using Priestley-Taylor:
```jldoctest; output=false
Tair, pressure, Rn = 30.0,100.0,500.0
ET_pot, LE_pot = potential_ET(PriestleyTaylor(), Tair, pressure, Rn)    
≈(ET_pot, 0.000204; rtol = 1e-2)
# output
true
``` 

Calculate potential ET for a surface with known Gs (0.5 mol m-2 s-1) and Ga (0.1 m s-1)
using Penman-Monteith:
```jldoctest; output=false
Tair, pressure, Rn = 30.0,100.0,500.0
VPD, Ga_h, Gs_pot = 2.0, 0.1, 0.5
ET_pot, LE_pot = potential_ET(
  PenmanMonteith(), Tair,pressure,Rn,VPD, Ga_h; Gs_pot)    
# now cross-check with the inverted equation
Gs_ms, Gs_mol = surface_conductance(
  InversePenmanMonteith(), Tair,pressure,VPD,LE_pot,Rn,Ga_h)
Gs_mol ≈ Gs_pot
# output
true
``` 

DataFrame variant with explicitly replacing missings using `coalesce.`:
```jldoctest; output=false
using DataFrames
df = DataFrame(
  Tair = 20.0:1.0:30.0,pressure = 100.0, Rn = 500.0, G = 105.0, VPD = 2.0, 
  Ga_h = 0.1) 
allowmissing!(df, Cols(:G)); df.G[1] = missing
#
# need to provide G explicitly
df_ET = potential_ET!(copy(df), PriestleyTaylor(); G = df.G)    
ismissing(df_ET.ET_pot[1])
#
# use coalesce to replace missing values by zero
df_ET = potential_ET!(
  copy(df), PriestleyTaylor(); G = coalesce.(df.G, zero(df.G)))    
!ismissing(df_ET.ET_pot[1])
# output
true
``` 
"""
function potential_ET(approach::PriestleyTaylor, Tair, pressure, Rn;
  G=zero(Tair),S=zero(Tair), kwargs...)
  #
  potential_ET(approach, Tair, pressure, Rn, G, S; kwargs...)
end
function potential_ET(::PriestleyTaylor, Tair, pressure, Rn, G, S;
  alpha=1.26,
  Esat_formula=Sonntag1990(),
  constants=BigleafConstants())
  #
  gamma  = psychrometric_constant(Tair,pressure;constants)
  Delta  = Esat_from_Tair_deriv(Tair; Esat_formula, constants)
  LE_pot = (alpha * Delta * (Rn - G - S)) / (Delta + gamma)
  ET_pot = LE_to_ET(LE_pot,Tair)
  (ET_pot = ET_pot, LE_pot = LE_pot)
end
function potential_ET!(df::AbstractDataFrame, approach::PriestleyTaylor; 
  G=0.0,S=0.0, kwargs...) 
  # G or S might be vectors, hence cannot use ByRow
  ft = (args...) -> potential_ET.(Ref(approach), args...,G,S; kwargs...)
  transform!(df, Cols(:Tair, :pressure, :Rn) => ft => AsTable)
end

function potential_ET(approach::PenmanMonteith, Tair, pressure, Rn, VPD, Ga_h;
  G=zero(Tair),S=zero(Tair), kwargs...)
  #
  potential_ET(approach, Tair, pressure, Rn, VPD, Ga_h, G, S; kwargs...)
end
function potential_ET(::PenmanMonteith, Tair, pressure, Rn, VPD, Ga_h, G, S;
  Gs_pot=0.6,
  Esat_formula=Sonntag1990(),
  constants=BigleafConstants()
  )
  gamma  = psychrometric_constant(Tair,pressure;constants)
  Delta  = Esat_from_Tair_deriv(Tair; Esat_formula = Esat_formula,constants)
  Gs_pot = mol_to_ms(Gs_pot,Tair,pressure;constants)
  rho    = air_density(Tair,pressure;constants)
  LE_pot = (Delta * (Rn - G - S) + rho * oftype(VPD,constants.cp) * VPD * Ga_h) / 
    (Delta + gamma * (1 + Ga_h / Gs_pot))
  ET_pot = LE_to_ET(LE_pot,Tair)
  (ET_pot = ET_pot, LE_pot = LE_pot)
end
function potential_ET!(df::AbstractDataFrame, approach::PenmanMonteith; 
  G=0.0,S=0.0, kwargs...
  ) 
  # G or S might be vectors, hence cannot use ByRow
  ft = (args...) -> potential_ET.(Ref(approach), args..., G, S; kwargs...)
  transform!(df, Cols(:Tair, :pressure, :Rn, :VPD, :Ga_h) => ft => AsTable)
end


# function potential_ET(df, approach::PriestleyTaylor; 
#   G=missing,S=missing, infoGS=true, kwargs...) 
#   #
#   dfGS = get_df_GS(df, G,S; infoGS) 
#   f(args...) = potential_ET(args..., approach; kwargs...)
#   select(hcat(select(df,:Tair, :pressure, :Rn), dfGS; copycols = false),
#     All() => ByRow(f) => AsTable
#   )
# end
# function potential_ET(df, approach::PenmanMonteith; 
#   G=missing,S=missing, infoGS=true, kwargs...) 
#   #
#   dfGS = get_df_GS(df, G,S; infoGS) 
#   function f(args...) 
#     potential_ET(args..., approach; kwargs...)
#   end
#   select(hcat(select(df,:Tair, :pressure, :Rn, :VPD, :Ga), dfGS; copycols = false),
#     All() => ByRow(f) => AsTable
#   )
# end

# function get_df_GS(df, G, S; infoGS=true)
#   nout = nrow(df)
#   G_ = if ismissing(G)
#     infoGS && @info("Ground heat flux G is not provided and set to 0.")
#     Zeros(nout)
#   else 
#     G
#   end
#   S_ = if ismissing(S)
#     infoGS && @info("Storage heat flux S is not provided and set to 0.")
#     Zeros(nout)
#   else
#     S
#   end
#   (G  = G_, S = S_)
# end

# function fill_GS_missings(df::DataFrame, G,S, missing_G_as_NA, missing_S_as_NA; infoGS=true)
#   nout = nrow(df)
#   G_ = ifelse(
#     ismissing(G),
#     (infoGS && @info("Ground heat flux G is not provided and set to 0."); Zeros(nout)),
#     @. ifelse(missing_G_as_NA, G, coalesce(G, 0.0))
#     )
#   S_ = ifelse(
#   ismissing(S),
#   (infoGS && @info("Storage heat flux S is not provided and set to 0."); Zeros(nout)),
#   @. ifelse(missing_S_as_NA, G, coalesce(S, 0.0))
#   )
#   DataFrame(G  = G_, S = S_)
# end

# function fill_vec(G; is_replace_missing = true, fillvalue = zero(G))
#   @. ifelse(is_replace_missing, coalesce(G, fillvalue), G)
# end

"""
    equilibrium_imposed_ET(Tair,pressure,VPD,Gs, Rn; ...)
    equilibrium_imposed_ET!(df; ...)

Evapotranspiration (ET) split up into imposed ET and equilibrium ET.

# Arguments
- `Tair`      : Air temperature (deg C)
- `pressure`  : Atmospheric pressure (kPa)
- `VPD`       : Air vapor pressure deficit (kPa)
- `Gs`        : surface conductance to water vapor (m s-1)
- `Rn`        : Net radiation (W m-2)
optional : 
- `G=0`       : Ground heat flux (W m-2)
- `S=0`       : Sum of all storage fluxes (W m-2)
- `Esat_formula=Sonntag1990()`: formula used in [`Esat_from_Tair`](@ref)
- `constants=`[`BigleafConstants`](@ref)`()`: physical constants (cp, eps)
                 
# Details
Total evapotranspiration can be written in the form (Jarvis & McNaughton 6):

``ET = \\Omega \\mathit{ET}_{eq} + (1 - \\Omega) \\mathit{ET}_{imp}``

where ``\\Omega`` is the decoupling coefficient as calculated from
[`decoupling`](@ref). `ET_eq` is the equilibrium evapotranspiration i.e.,
the ET rate that would occur under uncoupled conditions, where the budget
is dominated by radiation (when Ga -> 0):

``ET_{eq} = (\\Delta \\, (R_n - G - S)) / (\\lambda \\, (\\Delta + \\gamma))``

where ``\\Delta`` is the slope of the saturation vapor pressur(kPa K-1),
``\\lambda`` is the latent heat of vaporization (J kg-1), and ``\\gamma``
is the psychrometric constant (kPa K-1).
`ET_imp` is the imposed evapotranspiration rate, the ET rate
that would occur under fully coupled conditions (when Ga -> inf):

``ET_{imp} = (\\rho \\, c_p \\, \\mathit{VPD} ~ G_s) / (\\lambda \\, \\gamma)``

where ``\\rho`` is the air density (kg m-3).
 
# Value
A `NamedTuple` or `DataFrame` with the following columns:
- `ET_eq`: Equilibrium ET (kg m-2 s-1)
- `ET_imp`: Imposed ET (kg m-2 s-1)
- `LE_eq`: Equilibrium LE (W m-2)
- `LE_imp`: Imposed LE (W m-2)      

# References
- Jarvis, P_G., McNaughton, K_G., 1986: Stomatal control of transpiration:
  scaling up from leaf to region. Advances in Ecological Rese1-49.
- Monteith, J_L., Unsworth, M_H., 2008: Principles of ironmPhysics.
  3rd edition. Academic Press, London. 
            
# Examples
```jldoctest; output = false
Tair,pressure,Rn, VPD, Gs = 20.0,100.0,50.0, 0.5, 0.01
ET_eq, ET_imp, LE_eq, LE_imp = equilibrium_imposed_ET(Tair,pressure,VPD,Gs, Rn)    
≈(ET_eq, 1.399424e-05; rtol = 1e-5)
# output
true
``` 
"""
function equilibrium_imposed_ET(Tair,pressure,VPD,Gs, Rn;
  G=zero(Tair),S=zero(Tair), kwargs...)
  equilibrium_imposed_ET(Tair,pressure,VPD,Gs, Rn, G, S; kwargs...)
end

function equilibrium_imposed_ET(Tair,pressure,VPD,Gs, Rn, G, S;
  Esat_formula=Sonntag1990(),
  constants=BigleafConstants())
  # 
  rho    = air_density(Tair, pressure; constants)
  gamma  = psychrometric_constant(Tair, pressure; constants)
  Delta  = Esat_from_Tair_deriv(Tair; Esat_formula = Esat_formula, constants)
  LE_eq  = (Delta * (Rn - G - S)) / (gamma + Delta)
  LE_imp = (rho * constants.cp * Gs * VPD) / gamma
  #
  ET_imp = LE_to_ET(LE_imp,Tair)
  ET_eq  = LE_to_ET(LE_eq,Tair)
  (;ET_eq, ET_imp, LE_eq, LE_imp)
end
function equilibrium_imposed_ET!(df; G=0.0, S=0.0, kwargs...) 
  f(args...) = equilibrium_imposed_ET.(args..., G, S; kwargs...)
  transform!(df, [:Tair, :pressure, :VPD, :Gs, :Rn] => f => AsTable)
end

"""
    TODO; implement decoupling.

This stub is there to satisfy links im Help-pages.
"""
function decoupling()
  error("not yet implemented.")
end


