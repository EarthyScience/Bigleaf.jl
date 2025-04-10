abstract type SurfaceConductanceMethod end
struct FluxGradient <: SurfaceConductanceMethod end
struct InversePenmanMonteith <: SurfaceConductanceMethod end


"""
      surface_conductance(::FluxGradient,   Tair,pressure,VPD,LE; constants)
      surface_conductance(::InversePenmanMonteith, Tair,pressure,VPD,LE,Rn,Ga_h; 
            G = 0.0, S = 0.0, Esat_formula=Sonntag1990(), constants)

      surface_conductance!(df, method::FluxGradient; kwargs...)
      surface_conductance!(df, method::InversePenmanMonteith; G = 0.0, S = 0.0, kwargs...)

Calculate surface conductance to water vapor from the inverted Penman-Monteith equation 
or from a simple flux-gradient approach.

# Arguments
- `Tair`      : Air temperature (deg C)
- `pressure`  : Atmospheric pressure (kPa)
- `Rn`        : Net radiation (W m-2)
- `df`        : DataFrame with above variables
- `constants=`[`BigleafConstants`](@ref)`()`: Dictionary with physical constants

additional for InversePenmanMonteith
- `LE`        : Latent heat flux (W m-2)
- `VPD`       : Vapor pressure deficit (kPa)
- `Ga_h`      : Aerodynamic conductance towater vapor (m s-1), assumed equal to that of heat
- `G=0.0`     : Ground heat flux (W m-2), defaults to zero
- `S=0.0`     : Sum of all storage fluxes (W m-2), defaults to zero
- `Esat_formula=Sonntag1990()`: formula used in [`Esat_from_Tair`](@ref)

# Details
For `InversePenmanMonteith()`, surface conductance (Gs) in m s-1 
is calculated from the inverted Penman-Monteith equation:

``G_s = ( LE \\, G_a \\, \\gamma ) / ( \\Delta \\, A + \\rho \\, c_p \\, G_a \\, VPD - 
LE \\, (\\Delta + \\gamma ))``

Where ``\\gamma`` is the psychrometric constant (kPa K-1), ``\\Delta`` is the slope of the 
saturation vapor pressure curve (kPa K-1), and ``\\rho`` is air density (kg m-3).
Available energy (A) is defined as ``A = R_n - G - S``. 

Ground heat flux and total storage flux can be provided as scalars or vectors of
the length of the DataFrame in the DataFrame variant.
While the bigleaf R package by default converts any missings in `G` and `S` to 0,
in `Bigleaf.jl` the caller must take care, e.g. by using `G = coalesce(myGvector, 0.0)`.

For `FluxGradient()`, Gs (in mol m-2 s-1) is calculated from VPD and ET only:

   ``Gs = ET/p \\, VPD``

where ET is in mol m-2 s-1 and p is pressure. 
Note that this formulation assumes fully coupled conditions (i.e. Ga = inf).
This formulation is equivalent to the inverted form of Eq.6 in McNaughton & Black 1973:

   ``Gs = LE \\, \\gamma / (\\rho \\, c_p \\, VPD)``

which gives Gs in m s-1. Note that Gs > Gc (canopy conductance) under conditions 
when a significant fraction of ET comes from interception or soil evaporation. 

If `pressure` is not available, it can be approximated by elevation using the 
function [`pressure_from_elevation`](@ref)
 
# Value
NamedTuple with entries: 
 - Gs_ms: Surface conductance in m s-1
 - Gs_mol: Surface conductance in mol m-2 s-1

# References
- Monteith, J., 1965: Evaporation and environment. In Fogg, G. E. (Ed.),
  The state and movement of water in living organisms (pp.205-234).
  19th Symp. Soc. Exp. Biol., Cambridge University Press, Cambridge
- McNaughton, K_G., Black, T_A., 1973: A study of evapotranspiration
  from a Douglas Fir forest using the energy balance approach. Water
  Resources Research 9, 1579-1590.

# Examples
```jldoctest; output = false
Tair,pressure,VPD,LE,Rn,Ga_h,G = (14.8, 97.7, 1.08, 183.0, 778.0, 0.116, 15.6)
Gs = surface_conductance(InversePenmanMonteith(), Tair,pressure,VPD,LE,Rn,Ga_h;G)
isapprox(Gs.Gs_mol, 0.28, atol=0.1)
# output
true
``` 
"""
function surface_conductance(::FluxGradient, Tair,pressure,VPD,LE; constants=BigleafConstants())
    Gs_mol = (LE_to_ET(LE,Tair)/oftype(pressure,constants.Mw)) * pressure / VPD
    Gs_ms  = mol_to_ms(Gs_mol,Tair,pressure; constants)
    (;Gs_ms,Gs_mol)   
end     
function surface_conductance!(df::AbstractDataFrame, method::FluxGradient; kwargs...)
      fr = (args...) -> surface_conductance(method, args...; kwargs...) 
      transform!(df, SA[:Tair,:pressure,:VPD,:LE] => ByRow(fr) => AsTable)
end

function surface_conductance(::InversePenmanMonteith, Tair,pressure,VPD,LE,Rn,Ga_h; 
      G = 0.0, S = 0.0,
      #Ga="Ga_h",missing_G_as_NA=false,missing_S_as_NA=false,
      Esat_formula=Sonntag1990(), constants=BigleafConstants()
      )
      Delta = Esat_from_Tair_deriv(Tair; Esat_formula, constants)
      gamma = psychrometric_constant(Tair, pressure; constants)
      rho   = air_density(Tair, pressure; constants)
      Gs_ms  = ( LE * Ga_h * gamma ) / 
      (Delta * (Rn-G-S) + rho * oftype(Tair,constants.cp) * Ga_h * VPD - LE * (Delta + gamma))
      Gs_mol = ms_to_mol(Gs_ms,Tair,pressure; constants)
      (;Gs_ms,Gs_mol)    
end
function surface_conductance!(df::AbstractDataFrame, method::InversePenmanMonteith; 
      G = 0.0, S = 0.0, kwargs...)
      # need to convert G and S to positional arguments for proper broadcast
      fpos = (G, S, args...) -> surface_conductance(method, args...; G,S, kwargs...) 
      ft = (args...) -> fpos.(G, S, args...)
      transform!(df, SA[:Tair,:pressure,:VPD,:LE,:Rn,:Ga_h] => ft => AsTable)
end
    
