"""
    Reynolds_Number(Tair,pressure,ustar,z0m; constants)

calculates the Roughness Reynolds Number.

# Arguments
- `Tair`      : Air temperature (deg C)
- `pressure`  : Atmospheric pressure (kPa)
- `ustar`     : Friction velocity (m s-1)
- `z0m`       : Roughness length (m)
optional
- `constants=`[`BigleafConstants`](@ref)`()`
                 
# Details
The Roughness Reynolds Number is calculated as in Massman 1999a:     
``Re = z0m * ustar / v``, where `v` is the kinematic viscosity (m2 s-1).
         
# Value
Roughness Reynolds Number (-)

# References
Massman, W_J., 1999a: A model study of kB H- 1 for vegetated surfaces using
'localized near-field' Lagrangian theory. Journal of Hydrology 223, 27-43.

```jldoctest; output = false
Tair,pressure,ustar,z0m = 25.0,100.0,0.5,0.5
R = Reynolds_Number(Tair,pressure,ustar,z0m)                             
≈(R, 15870, rtol=1e-3) 
# output
true
``` 
"""
function Reynolds_Number(Tair,pressure,ustar,z0m;constants=BigleafConstants())
  v  = kinematic_viscosity(Tair,pressure;constants)
  Re = z0m*ustar/v
end


abstract type RoughnessMethod end
struct Roughness_wind_profile <: RoughnessMethod end
struct RoughnessCanopyHeight <: RoughnessMethod end
struct RoughnessCanopyHeightLAI <: RoughnessMethod end

"""
    roughness_parameters(::RoughnessCanopyHeight()    , zh; frac_d=0.7, frac_z0m=0.1)
    roughness_parameters(::RoughnessCanopyHeightLAI(), zh, LAI; cd=0.2, hs=0.01)
    roughness_parameters(::Roughness_wind_profile()     , ustar, wind, psi_m; 
      zh, zr, d = 0.7*zh, constants)

    roughness_parameters(::Roughness_wind_profile(), ustar, wind, Tair, pressure, H; 
      zh, zr, d = 0.7*zh, stab_formulation=Dyer1970(), constants)
    roughness_parameters(::Roughness_wind_profile(), df::DFTable; ...)
    

A approximations of the two roughness parameters displacement height (d)
and roughness length for momentum (z0m).
             
# Arguments              
- `zh`        : Vegetation height (m)          
- `constants=`[`BigleafConstants`](@ref)`()`: 

By canopy height:
- `frac_d`    : Fraction of displacement height on canopy height (-)
- `frac_z0m`  : Fraction of roughness length on canopy height (-)

By canopy height and LAI
- `LAI`       : Leaf area index (-) 
- `cd`        : Mean drag coefficient for individual leaves. Defaults to 0.2. 
- `hs`        : roughness length of the soil surface (m). 

By wind profile
- `wind`      : Wind speed at height zr (m s-1)
- `ustar`     : Friction velocity (m s-1)
- `zr`        : Instrument (reference) height (m)
- `d`         : Zero-plane displacement height (-)
- `psi_m`     : value of the stability function for heat, see [`stability_correction`](@ref)

Another variant estimates of `psi_m` by [`stability_correction`](@ref), which
requires further input arguments.
For convenience, these arguments can be provided using a DataFrame, however, the result
then is not type stable.

# Details
The two main roughness parameters, the displacement height (d)
and the roughness length for momentum (z0m) can be estimated from simple
empirical relationships with canopy height (zh). If `method = RoughnessCanopyHeight()`,
the following formulas are used:  

``d = frac_d * zh``

``z0m = frac_{z0m} * zh``

where ``frac_d`` defaults to 0.7 and ``frac_{z0m}`` to 0.1.

Alternatively, d and z0m can be estimated from both canopy height and LAI
(If `method = RoughnessCanopyHeightLAI()`).
Based on data from Shaw & Pereira 1982, Choudhury & Monteith 1988 proposed 
the following semi-empirical relations:

``X = cd * LAI``

``d = 1.1 * zh * ln(1 + X^{1/4})`` 

``z0m = hs + 0.3 * zh * X^{1/2}``   for ``0 <= X <= 0.2``

``z0m = hs * zh * (1 - d/zh)``   for ``0.2 < X`` 

If `method = Roughness_wind_profile()`, z0m is estimated by solving
the wind speed profile for z0m:

``z0m = median((zr - d) * exp(-k*wind / ustar - psi_m)``
        
By default, d in this equation is fixed to 0.7*zh, but can be set to any
other value.   

# Value
a NamedTuple with the following components:
- `d`: Zero-plane displacement height (m)
- `z0m`: Roughness length for momentum (m)
- `z0m_se`: Only if `method = wind_profile`: Standard Error of the median for z0m (m)

# References
- Choudhury, B. J., Monteith J_L., 1988: A four-layer model for the heat
  budget of homogeneous land surfaces. Q. J. R. Meteorol. Soc. 114, 373-398.
- Shaw, R. H., Pereira, A., 1982: Aerodynamic roughness of a plant canopy: 
  a numerical experiment. Agricultural Meteorology, 26, 51-65.
   
# See also
[`wind_profile`](@ref)
    
```jldoctest; output = false
using DataFrames
# estimate d and z0m from canopy height for a dense (LAI=5) and open (LAI=2) canopy
zh = 25.0
roughness_parameters(RoughnessCanopyHeightLAI(),zh,5)
roughness_parameters(RoughnessCanopyHeightLAI(),zh,2)   
   
# fix d to 0.7*zh and estimate z0m from the wind profile
df = DataFrame(
  Tair=[25.0,25.0,25.0],pressure=100.0,wind=[3.0,4.0,5.0],
  ustar=[0.5,0.6,0.65],H=200.0)
roughness_parameters(Roughness_wind_profile(),df;zh,zr=40,d=0.7*zh)

# assume d = 0.8*zh
rp = roughness_parameters(Roughness_wind_profile(),df;zh,zr=40,d=0.8*zh)
≈(rp.z0m, 0.55, rtol=0.1)
# output
true
```
"""
function roughness_parameters(::RoughnessCanopyHeight, zh::FT;
    frac_d=FT(0.7), frac_z0m=FT(0.1)) where FT
  d      = frac_d*zh
  z0m    = frac_z0m*zh
  z0m_se = missing
  (;d, z0m, z0m_se)
end

function roughness_parameters(::RoughnessCanopyHeightLAI, zh::FT, LAI;
  cd=FT(0.2), hs=FT(0.01)) where FT
  X = cd * LAI
  d = FT(1.1) * zh * log(1 + X^FT(1/4))
  z0m = ifelse(0 <= X <= FT(0.2), hs + FT(0.3) * X^(1/2), FT(0.3) * zh * (1 - d/zh))
  z0m_se = missing
  (;d, z0m, z0m_se)
end

function roughness_parameters(::Roughness_wind_profile, ustar::AbstractVector{FT}, wind, psi_m;
  zh, zr, d = 0.7*zh, constants=BigleafConstants()
  ) where FT
  FT == Missing && return((d = d, z0m = missing, z0m_se = missing))
  z0m_all = allowmissing(@. (nonmissingtype(FT)(zr) - nonmissingtype(FT)(d)) * 
    exp(-nonmissingtype(FT)(constants.k)*wind / ustar - psi_m))
  #z0m_all[(z0m_all .> zh)] .= missing # problems with missings
  replace!(x -> !ismissing(x) && x > zh ? missing : x, z0m_all)
  nval = sum(.!ismissing.(z0m_all))
  z0m    = median(skipmissing(z0m_all))
  z0m_se0 = (std(skipmissing(z0m_all)) / sqrt(nval))
  # http://influentialpoints.com/Training/standard_error_of_median.htm
  z0m_se = oftype(z0m_se0,constants.se_median) * z0m_se0
  (;d, z0m, z0m_se)
end

function roughness_parameters(method::Roughness_wind_profile, 
  ustar::AbstractVector, wind, Tair, pressure, H; zh, zr, d = 0.7*zh,
  stab_formulation=Dyer1970(), constants=BigleafConstants(), kwargs...
  )
  # psi_m = Tables.Columns(stability_correction.(
  #   zr, d, Tair, pressure, ustar, H; stab_formulation, constants)).psi_m
  psi_m = getindex.(stability_correction.(
    zr, d, Tair, pressure, ustar, H; stab_formulation, constants), :psi_m)
  # https://discourse.julialang.org/t/eltype-changed-during-broadcast-with-missing/71312?u=earthyscience
  res_eltype = Union{eltype(ustar),eltype(wind),eltype(Tair),eltype(pressure),eltype(H)}
  psi_m = convert(Vector{res_eltype}, psi_m)::Vector{res_eltype}    
  rp = roughness_parameters(method, ustar, wind, psi_m; zh, zr, d, constants, kwargs...)
end

function roughness_parameters(method::Roughness_wind_profile, df::DFTable; 
  psi_m = nothing, stab_formulation=Dyer1970(), kwargs...
  )
  !isnothing(psi_m) && return(roughness_parameters(
    method, df.ustar, df.wind, psi_m; kwargs...))
  if stab_formulation isa NoStabilityCorrection
    psi_m = 0.0
    roughness_parameters(method, df.ustar, df.wind, psi_m; kwargs...)
  else
    roughness_parameters(method, df.ustar, df.wind, df.Tair, df.pressure, df.H; 
      stab_formulation, kwargs...)
  end
end

 
"""                                                                                                                       
    wind_profile(z::Number, ustar, d, z0m, psi_m = zero(z); constants)
    wind_profile(z, df::AbstractDataFrame, d, z0m, psi_m::AbstractVector; constants)
    wind_profile(z, df::DFTable, d, z0m; psi_m = nothing,
      stab_formulation = Dyer1970(), constants)

Wind speed at a given height above the canopy estimated from single-level
measurements of wind speed.

# Arguments
- `z`         : Height above ground for which wind speed is calculated.
- `ustar`     : Friction velocity (m s-1)
- `d`         : Zero-plane displacement height (-)
- `z0m`       : Roughness length (m)
- `constants=`[`BigleafConstants`](@ref)`()`
For DataFrame variant with supplying stability_parameter
- `df`:         : DataFrame with columns 
  - `ustar`     : Friction velocity (m s-1)
- `psi_m`     : value of the stability function for heat, see [`stability_correction`](@ref)
  Pass `psi_m = 0.0` to neglect stability correction.
For DataFrame varinat where psi_m is to be estimated
- additional columns in `df`:
  - `Tair`, `pressure`, `H` : see [`stability_correction`](@ref)

# Details
The underlying assumption is the existence of a logarithmic wind profile
above the height d + z0m (the height at which wind speed mathematically reaches zero
according to the Monin-Obhukov similarity theory).
In this case, the wind speed at a given height z is given by:

``u(z) = (ustar/k) * (ln((z - d) / z0m) - \\psi_m``

The roughness parameters zero-plane displacement height (d) and roughness length (z0m)
can be approximated from [`roughness_parameters`](@ref). 
``\\psi_m`` is the stability correction. Set it to zero (not recommended) to neglect
statbility correction. By default it is estimated from wind profile using
[`stability_correction`](@ref)
                                                            
# Note
Note that this equation is only valid for z >= d + z0m, and it is not 
meaningful to calculate values closely above d + z0m. All values in `heights`
smaller than d + z0m will return 0. 

                              
# Value
wind speed at given height `z`. 
        
# References
- Monteith, J_L., Unsworth, M_H., 2008: Principles of Environmental Physics.
  3rd edition. Academic Press, London. 
- Newman, J_F., Klein, P_M., 2014: The impacts of atmospheric stability on
  the accuracy of wind speed extrapolation methods. Resources 3, 81-105.
        
# See also
[`roughness_parameters`](@ref)

```jldoctest; output = false
using DataFrames
heights = 18:2:40  # heights above ground for which to calculate wind speed
df = DataFrame(
  Tair=25.0,pressure=100.0,wind=[3.0,4,5],ustar=[0.5,0.6,0.65],H=[200.0,230,250]) 
zr=40;zh=25;d=16
# z0m and MOL are independent of height, compute before
MOL = Monin_Obukhov_length.(df.Tair, df.pressure, df.ustar, df.H)
z0m = roughness_parameters(
  Roughness_wind_profile(), df.ustar, df.wind, df.Tair, df.pressure, df.H; zh, zr).z0m 
ws = map(heights) do z
  wind_profile(z,df,d,z0m; MOL)
end
using Plots # plot wind profiles for the three rows in df
plot(first.(ws), heights, ylab = "height (m)", xlab = "wind speed (m/s)", legend=:topleft)
plot!(getindex.(ws, 2), heights)
plot!(getindex.(ws, 3), heights)
nothing
# output
``` 
"""  
function wind_profile(
  z::Number, ustar::Union{Missing,FT}, d, z0m, psi_m; constants=BigleafConstants()
  ) where FT<:Number
  wind_heights = max(FT(0),(ustar / FT(constants.k)) * (log(max(FT(0),(FT(z) -FT(d))) / FT(z0m)) - psi_m))
  wind_heights
end

function wind_profile(z::Number, ustar::Union{Missing,Number}, d, z0m, Tair,pressure,H,
  stab_formulation=Dyer1970(), constants=BigleafConstants())
  psi_m = stability_correction(
    z,d, Tair,pressure,ustar,H; stab_formulation, constants).psi_m
  wind_profile(z, ustar, d, z0m, psi_m)
end

function wind_profile(z::Number, df::DFTable, d, z0m; psi_m = nothing, MOL = nothing,
  stab_formulation = Dyer1970(), constants =BigleafConstants()
  )
  if isnothing(psi_m)
    if isnothing(MOL) 
      psi_m = stability_correction(df; z, d, stab_formulation, constants).psi_m
    else
      zeta = stability_parameter.(z,d,MOL)
      psi_m = getindex.(stability_correction.(zeta), :psi_m)
    end
  end
  #wind_profile(df, z, d, z0m, psi_m; constants)
  windz = wind_profile.(z, df.ustar, d, z0m, psi_m; constants)
end

