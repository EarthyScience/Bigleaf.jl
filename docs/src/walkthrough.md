This vignette is a short introduction to the functionalities of the `Bigleaf.jl` package. 
It is directed to first-time package users who are familiar with the basic concepts of Julia. 
After presenting the use of several key functions of the package, 
some useful hints and guidelines are given at the end of the vignette.


# Package scope and important conceptual considerations

`Bigleaf.jl` calculates physical and physiological ecosystem properties from eddy 
covariance data. Examples for such properties are aerodynamic and surface conductance, 
surface conditions(e.g. temperature, VPD), wind profile, roughness parameters, 
vegetation-atmosphere decoupling, potential evapotranspiration, (intrinsic) water-use 
efficiency, stomatal sensitivity to VPD, or intercellular CO2 concentration.  
All calculations in the `Bigleaf.jl` package assume that the ecosystem behaves like a  
"big-leaf", i.e. a single, homogeneous plane which acts as the only source and sink of the 
measured fluxes. This assumption comes with the advantages that calculations are simplified 
considerably and that (in most cases) little ancillary information on the EC sites is 
required. It is important to keep in mind that these simplifications go hand in hand with 
critical limitations. All derived variables are bulk ecosystem characteristics and have to 
be interpreted as such. It is for example not possible to infer within-canopy variations 
of a certain property.

Please also keep in mind that the `Bigleaf.jl` package does NOT provide formulations for 
bottom-up modelling. The principle applied here is to use an inversion approach in which 
ecosystem properties are inferred top-down from the measured fluxes. Such an inversion can, 
in principle, be also be conducted with more complex models (e.g. sun/shade or canopy/soil 
models), but keep in mind that these approaches also require that the additional, 
site-specific parameters are adequately well known. 

The use of more detailed models is not within the scope of the `Bigleaf.jl` package, but 
it is preferable to use such approaches when important assumptions of the "big-leaf" 
approach are not met. This is the case in particular when the ecosystem is sparsely covered 
with vegetation (low LAI, e.g. sparse crops, some savanna systems). 


# Preparing the data

In this tutorial, we will work with a dataset from the eddy covariance site Tharandt 
(DE-Tha), a spruce forest in Eastern Germany. The DataFrame `DE_Tha_Jun_2014` is downloaded 
from the `bigleaf` 
[R package](https://bitbucket.org/juergenknauer/Bigleaf/) repository and contains 
half-hourly data of meteorological and flux measurements made in June 2014. For loading the 
RData into Julia, see the 
[source](https://github.com/earthyscience/Bigleaf.jl/blob/main/docs/src/walkthrough.md?plain=1#L26) 
of this file. We give the data.frame a shorter name here and create a timestamp.

```@example doc
using Bigleaf
using DataFrames
```
```@setup doc
using Latexify
using DataDeps, Suppressor
using RData
import CodecBzip2, CodecXz
#@suppress_err # error in github-actions: GitHubActionsLogger has no field stream
register(DataDep(
    "DE_Tha_Jun_2014.rda",
    "downloading example dataset DE_Tha_Jun_2014 from bitbucket.org/juergenknauer/bigleaf",
    "https://bitbucket.org/juergenknauer/bigleaf/raw/0ebe11626b4409305951e8add9f6436703c82584/data/DE_Tha_Jun_2014.rda",
    "395f02e1a1a2d175ac7499c200d9d48b1cb58ff4755dfd2d7fe96fd18258d73c"
))
#println(datadep"DE_Tha_Jun_2014.rda")
ENV["DATADEPS_ALWAYS_ACCEPT"]="true" # avoid question to download
DE_Tha_Jun_2014 = first(values(load(joinpath(
  datadep"DE_Tha_Jun_2014.rda/DE_Tha_Jun_2014.rda"))))
nothing
```
```@example doc
tha = DE_Tha_Jun_2014
set_datetime_ydh!(tha)
# mdtable(select(describe(tha), :variable, :eltype, :min, :max), latex=false) # hide
nothing # hide
```

And the first six rows of tha:
```@example doc
mdtable(tha[1:6,:],latex=false) # hide
```

More information on the data (e.g. meaning of column names and units) can be found at the 
[bigleaf R package](https://bitbucket.org/juergenknauer/bigleaf/src/master/man/DE_Tha_Jun_2014.Rd). 
For more information on the site see e.g. Grünwald & Bernhofer 2007.
In addition, we will need some ancillary data for this site throughout this tutorial. To ensure consistency, we define them here at the beginning:

```@example doc
thal = (
     LAI = 7.6,   # leaf area index
     zh  = 26.5,  # average vegetation height (m)
     zr  = 42,    # sensor height (m)
     Dl  = 0.01,  # leaf characteristic dimension (m)
)
nothing # hide
```

# General guidelines on package usage

There are a few general guidelines that are important to consider when using the `Bigleaf.jl` package. 


## Units

It is imperative that variables are provided in the right units, as the plausibility of 
the input units is not checked in most cases. The required units of the input arguments 
can be found in the respective help file of the function. The good news is that units 
do not change across functions. For example, pressure is always required in kPa, 
and temperature always in °C.

## Function arguments

`Bigleaf.jl` usually provides functions in two flavours.
- providing all arguments separately as scalars and output being a single scalar
  or a NamedTuple
- providing a DataFrame as first argument with columns corresponding to the inputs and 
  output being the in-place modified DataFrame. Most keyword arguments
  accept both, vectors or scalars.
  The column names in the DataFrame should correspond to the argument names
  of the corresponding method with individual inputs.

We can demonstrate the usage with a simple example:

```@example doc
# explicit inputs
Tair, pressure, Rn, =  14.81, 97.71, 778.17 
potential_ET(PriestleyTaylor(), Tair, pressure, Rn)

# DataFrame
potential_ET!(copy(tha), PriestleyTaylor())

# DataFrame with a few columns overwritten by user values
potential_ET!(transform(tha, :Tair => x -> 25.0; renamecols=false), PriestleyTaylor())

# varying one input only: scalar form with dot-notation
Tair_vec =  10.0:1.0:20.0
DataFrame(potential_ET.(Ref(PriestleyTaylor()), Tair_vec, pressure, Rn))
nothing # hide
```

For functions operating only on vectors, e.g. [`roughness_parameters`](@ref) vectors
are provided with the individual inputs. Sometimes, an additional  non-mutating DataFrame 
variant is provided for convenience, however, the output value of this variant is 
not type-stable.

## Ground heat flux and storage fluxes

Many functions require the available energy ($A$), which is defined as ($A = R_n - G - S$, 
all in $\text{W m}^{-2}$), where $R_n$ is the net radiation, $G$ is the ground heat flux, 
and $S$ is the sum of all storage fluxes of the ecosystem 
(see e.g. Leuning et al. 2012 for an overview). For some sites, $G$ is not available, 
and for most sites, only a few components of $S$ are measured. 

In `Bigleaf.jl` it is not a problem if $G$ and/or $S$ are missing (other than the results 
might be (slightly) biased), but special options exist for the treatment of missing 
$S$ and $G$ values. 

Note that the default for G and S in the dataframe variant is missing (and assumed zero), 
even if those columns are
present in the DataFrame. You need to explicitly pass those columns with the optional
arguments: e.g. `potential_ET(df, PriestleyTaylor(); G = df.G)`

Note that in difference to the bigleaf R package missing entries in an input
vector are not relaced by zero by default. 
You need to explicitly use coalesce when specifying a ground heat flux
for which missings should be replaced by zero: `;G = coalesce(df.G, zero(df.G))`
 
# Function walkthrough #

## Data filtering

For most applications it is meaningful to filter your data. There are two main reasons 
why we want to filter our data before we start calculating ecosystem properties. 
The first one is to exclude datapoints that do not fulfill the requirements of the 
EC technique or that are of bad quality due to e.g. instrument failure or gap-filling 
with poor confidence. Note that the quality assessment of the EC data is not the purpose 
of the `Bigleaf.jl` package. This is done by other packages (e.g. `REddyProc`), 
which often provide quality control flags for the variables. These quality control 
flags are used here to filter out bad-quality datapoints.

A second reason for filtering our data is that some derived properties are only 
meaningful if certain meteorological conditions are met. For instance, if we are 
interested in properties related to plant gas exchange, it makes most sense to focus on 
time periods when plants are photosynthetically active 
(i.e. in the growing season and at daytime).

`Bigleaf.jl` provides methods that update (or create) the :valid column in 
a DataFrame. Records, i.e. rows, that contain non valid conditions are set to false.
If the valid column was false before, it stays at false.

### `setinvalid_qualityflag!`
One can check quality flags. By default (argument `setvalmissing = true`) this also
replaces the non-valid values in the data-columns by `missing`.
```@example doc
thaf = copy(tha)   # keep the original tha DataFrame
# if the :valid columns does not exist yet, it is created with all values true
setinvalid_qualityflag!(thaf; vars = ["LE", "NEE"])
sum(.!thaf.valid) # 7 records marked non-valid
sum(ismissing.(thaf.NEE)) # 7 NEE values set to missing
```
In the function call above, `vars` lists the variables that should be filtered with 
respect to their quality. Optional parameter `qc_suffix="_qc"` denotes the extension 
of the variable name that identifies the column as a quality control indicator of a given 
variable. The variables `LE` and `LE_qc`, for example, denote the variable itself 
(latent heat flux), and the quality of the variable `LE`, respectively. The optional 
argument `good_quality_threshold = 1.0` specifies the values of the quality column
below which the quality control to be considered as acceptable quality 
(i.e. to not be filtered). For example with default value 1, 
all `LE` values whose `LE_qc` variable is larger than 1 are set to `missing`. 
The variable `missing_qc_as_bad` is required to decide what to do in 
case of missing values in the quality control variable. By default this is (conservatively) 
set to `true`, i.e. all entries where the qc variable is missing is set invalid. 

### `setinvalid_range!`

We can  filter for meteorological conditions to be in acceptable ranges. 
For each variable to check we supply the valid minimum and valid maximum as a two-tuple
as the second component of a pair. If their is no limit towards small or
large values, supply `-Inf` or `Inf` as the minimum or maximum respectively.
```@example doc
setinvalid_range!(thaf, 
     :PPFD => (200, Inf), 
     :ustar => (0.2, Inf), 
     :LE =>(0, Inf), 
     :VPD => (0.01, Inf)
     )
sum(.!thaf.valid) # many more records marked invalid
minimum(skipmissing(thaf.PPFD)) >= 200 # values outsides range some set to missing
sum(ismissing.(thaf.PPFD))
```

About half of the data were filtered because radiation was not high enough (night-time). 
Another quarter was filtered because they showed negative LE values. 
However, most of them occurred during the night:
```@example doc
sum(ismissing.(thaf.PPFD)) / nrow(thaf) # 0.48
sum(.!ismissing.(thaf.PPFD) .&& ismissing.(thaf.LE)) / nrow(thaf) # 0.05
```

### `setinvalid_nongrowingseason!`

A third method filters periods outside the growing season:
```@example doc
setinvalid_nongrowingseason!(thaf, 0.4) 
sum(.!thaf.valid) # tha dataset is all within growing season - no additional invalids
```

This function implements a simple growing season filter based on daily smoothed GPP time 
series. 
Arguments  `tGPP` determines how high daily GPP has to be in relation to its peak value 
within the year. In this case, the value of 0.4 denotes that smoothed GPP has to be at 
least 40% of the 95th quantile. 
Argument `ws` controls the degree of smoothing in the timeseries 
and should be between 10-20 days. 
The purpose of which is to minimize the high variation of GPP between days,
Argument `min_int` is a parameter that avoids that data are switching from 
inside the growing season and out from one day to the next. 
It determines the minimum number of days that a season should have. 
The growing season filter is applicable to all sites, with one more more growing seasons,
but its advisable that site-specific parameter settings are used.

In this example, it does not really make sense to filter for growing season, 
since it uses only one month of data of which we know that vegetation is active at the site. 
The algorithm realizes that and does not mark any additional data as invalid.

### `setinvalid_afterprecip!`

As a last step we will filter for precipitation events. 
This is often meaningful for ecophysiological studies because data during and shortly 
after rainfall events do not contain much information on the physiological activity 
of the vegetation because they comprise significant fractions of evaporation from the 
soil and plant surfaces. The purpose of such a filter is mostly to minimize the fraction 
of soil and interception evaporation on the total water flux. This filter simply excludes 
periods following a precipitation event. A precipitation event, here, is defined as any time 
step with a recorded precipitation higher than `min_precip` (in mm per timestep). 
The function then filters all time periods following a precipitation event. 
The number of subsequent time periods excluded is controlled by the argument `precip_hours`. 
Here, we exclude rainfall events and the following 24 hours.
The timestamps in the DataFrame must be sorted in increasing order.

```@example doc
setinvalid_afterprecip!(thaf; min_precip=0.02, hours_after=24)
sum(.!thaf.valid) # some more invalids
```

In this walkthrough we use the data as filtered above:
```@example doc
thas = subset(thaf, :valid)
nrow(thas)
```

With first 6 rows:
```@example doc
mdtable(thas[1:6,:],latex=false) # hide
```

## Meteorological variables

The `Bigleaf.jl` package provides calculation routines for a number of meteorological variables, which are basic to the calculation of many other variables. A few examples on their usage are given below:

```@example doc
# Saturation vapor pressure (kPa) and slope of the saturation vapor pressure curve (kPa K-1)
Esat_slope(25.0)
```
```@example doc
# psychrometric constant (kPa K-1)
psychrometric_constant(25.0,100.0) # Tair, pressure
```
```@example doc
# air density (kg m-3)
air_density(25.0,100.0) # Tair, pressure
```
```@example doc
# dew point (degC)
dew_point(25.0,1.0) # Tair, VPD
```
```@example doc
# wetbulb temperature (degC)
wetbulb_temp(25.0, 100.0, 1.0) # Tair, pressure, VPD
```
```@example doc
# estimate atmospheric pressure from elevation (hypsometric equation)
pressure_from_elevation(500.0, 25.0) # elev, Tair
```

There are several formulations describing the empirical function `Esat(Tair)`.
The following figure compares them at absole scale and as difference to the 
#default method. The differences are small.

```@setup doc
#using DataFrames
#Tair = 0:0.25:12
##Tair = [10.0,20.0]
#eform_def = Sonntag1990()
#Esat_def = Esat_from_Tair.(Tair; Esat_formula = eform_def)
#eforms = (Sonntag1990(), Alduchov1996(), Allen1998())
#eform = eforms[2]
#string.(eforms)
#df = mapreduce(vcat, eforms) do eform 
#    Esat = Esat_from_Tair.(Tair; Esat_formula = eform)
#    local dff # make sure to not override previous results
#    dff = DataFrame(
#        Esat_formula = eform, Tair = Tair, 
#        Esat = Esat,
#        dEsat = Esat - Esat_def,
#        )
#end;
##using Chain
#using Pipe
#using Plots, StatsPlots
#dfw = @pipe df |> select(_, 1,2, :Esat) |> unstack(_, :Esat_formula, 3)
#dfws = @pipe df |> select(_, 1,2, :dEsat) |> unstack(_, :Esat_formula, 3)
#@df dfw plot(:Tair, cols(2:4), legend = :topleft, xlab="Tair (degC)", #ylab="Esat (kPa)")
#savefig("Esat_abs.svg")
#@df dfws plot(:Tair, cols(2:4), legend = :topleft, xlab="Tair (degC)", #ylab="Esat -ESat_Sonntag1990 (kPa)")
#savefig("fig/Esat_rel.svg")
```

![](fig/Esat_abs.svg)

![](fig/Esat_rel.svg)

## Aerodynamic conductance

An important metric for many calculations in the `Bigleaf.jl` package is the aerodynamic 
conductance ($G_a$) between the land surface and the measurement height. $G_a$ 
characterizes how efficiently mass and energy is transferred between the land surface 
and the atmosphere. $G_a$ consists of two parts: $G_{a_m}$, the aerodynamic conductance 
for momentum, and $G_b$, the canopy boundary layer (or quasi-laminar) conductance. 
$G_a$ can be defined as 

  $G_a = 1/(1/G_{a_m} + 1/G_b)$. 

In this tutorial we will focus on 
how to use the function [`aerodynamic_conductance!`](@ref). 
For further details on the equations, 
the reader is directed to the publication of the Bigleaf package (Knauer et al. 2018) and 
the references therein. A good overview is provided by e.g. Verma 1989.

  $G_a$ and in particular $G_b$ can be calculated with varying degrees of complexity. 
We start with the simplest version, in which $G_b$ is calculated empirically based on 
the friction velocity ($u_*$) according to Thom 1972:

```@example doc
aerodynamic_conductance!(thas);
thas[1:3, Cols(:datetime,Between(:zeta,:Ga_CO2))]
```

Note that by not providing additional arguments, the default values are taken.
We also do not need most of the arguments that can be provided to the function in this case 
(i.e. with `Gb_model=Thom1972()`). These are only required if we use a more complex 
formulation of $G_b$.
The output of the function is another DataFrame which contains separate columns for 
conductances and resistances of different scalars (momentum, heat, and $CO_2$ by default).

For comparison, we now calculate a second estimate of $G_a$, where the calculation of 
$G_b$ is more physically-based (Su et al. 2001), and which requires more input variables 
compared to the first version. In particular, we now need LAI, the leaf characteristic 
dimension ($D_l$, assumed to be 1cm here), and information on sensor and canopy height 
($z_r$ and $z_h$), as well as the displacement height (assumed to be 0.7*$z_h$):


```@example doc
aerodynamic_conductance!(thas;Gb_model=Su2001(),
     LAI=thal.zh, zh=thal.zh, d=0.7*thal.zh, zr=thal.zr, Dl=thal.Dl);
thas[1:3, Cols(:datetime,Between(:zeta,:Ga_CO2))]
```

We see that the values are different compared to the first, empirical estimate. 
This is because this formulation takes additional aerodynamically relevant properties 
(LAI, $D_l$) into account that were not considered by the simple empirical formulation.


## Boundary layer conductance for trace gases

By default, the function `aerodynamic_conductance` (calling `compute_Gb!`) returns the 
(quasi-laminar) canopy boundary layer ($G_{b}$) for heat and water vapor 
(which are assumed to be equal in the `Bigleaf.jl`), as well as for CO$_2$. 
Function `add_Gb` calculates $G_b$ for other trace gases, provided that the respective Schmidt 
number is known. 

```@example doc
compute_Gb!(thas, Thom1972()); # adds/modifies column Gb_h and Gb_CO2
add_Gb!(thas, :Gb_O2 => 0.84, :Gb_CH4 => 0.99); # adds Gb_O2 and Gb_CH4
select(first(thas,3), r"Gb_")
```

## Surface conductance

Knowledge of aerodynamic conductance $G_a$ 
allows us to calculate the bulk surface conductance ($G_s$) of the site 
(In this case by inverting the Penman-Monteith equation). Gs represents the combined 
conductance of the vegetation and the soil to water vapor transfer (and as such it is not 
a purely physiological quantity). Calculating $G_s$ in `Bigleaf.jl` is simple:

```@example doc
surface_conductance!(thas, InversePenmanMonteith());
thas[1:3,Cols(:datetime, r"Gs")]
```

The two columns only differ in the unit of $G_s$. 
One in m s$^{-1}$ and one in mol m$^{-2}$ s$^{-1}$. 
In this function we have ignored the ground heat flux ($G$) and the storage fluxes ($S$).
By default they are assumed zero.
In our example we do not have information on the storage fluxes, but we have measurements 
on the ground heat flux, which we should add to the function call:

```@example doc
surface_conductance!(thas, InversePenmanMonteith(); G=thas.G);
thas[1:3,Cols(:datetime, r"Gs")]
```


## Wind profile

The 'big-leaf' framework assumes that wind speed is zero at height d + $z_{0m}$ 
(where $z_{0m}$ is the roughness length for momentum) and then increases exponentially with 
height. The shape of the wind profile further depends on the stability conditions of the 
air above the canopy.
In `Bigleaf.jl`, a wind profile can be calculated assuming an exponential increase with 
height, which is affected by atmospheric stability. Here, we calculate wind speed at 
heights of 22-60m in steps of 2m. As expected, the gradient in wind speed is strongest 
close to the surface and weaker at greater heights:

```@example doc
using Statistics
wind_heights = 22:2:60.0
d = 0.7 * thal.zh
z0m = roughness_parameters(Roughness_wind_profile(), thas; zh=thal.zh, zr=thal.zr).z0m
wp = map(wind_heights) do z
  wind_profile(z, thas,d, z0m)
end;
nothing # hide
```
```@setup doc
wp_means = map(x -> mean(skipmissing(x)), wp)
wp_sd    = map(x -> std(skipmissing(x)), wp)
wr_mean = mean(skipmissing(thas.wind)) # measurements at reference height
wr_sd    = std(skipmissing(thas.wind))
using Plots # plot wind profiles for the three rows in df
plot(wp_means, wind_heights, ylab = "height (m)", xlab = "wind speed (m/s)", xerror=wp_sd, 
  label=nothing)
scatter!(wp_means, wind_heights, label = nothing)
```
```@example doc
scatter!([wr_mean], [thal.zr], xerror = [wr_sd], markerstrokecolor=:blue, #hide
  markerstrokewidth=2, label = nothing) # hide
```

Here, the points denote the mean wind speed and the bars denote the standard deviation
across time. The blue point/bar represent the values that were measured at zr = 42m. 
In this case we see that the wind speed as "back-calculated" from the wind profile agrees 
well with the actual measurements.


## Potential evapotranspiration

For many hydrological applications, it is relevant to get an estimate on the potential 
evapotranspiration (PET). At the moment, the `Bigleaf.jl` contains two formulations 
for the estimate of PET: the Priestley-Taylor equation, and the Penman-Monteith equation:

```@example doc
potential_ET!(thas, PriestleyTaylor(); G = thas.G)

# aerodynamic Ga_h and surface conductance Gs_mol must be computed before
potential_ET!(thas, PenmanMonteith();  G = thas.G, 
        Gs_pot=quantile(skipmissing(thas.Gs_mol),0.95))
thas[24:26, Cols(:datetime, :ET_pot, :LE_pot)]
```

In the second calculation it is important to provide an estimate of aerodynamic 
conductance, ``G_a``, and the potential surface conductance under optimal conditions, 
``G_{s pot}``. 
Here, we have approximated ``G_{s pot}`` with the ``95^{\text{th}}`` percentile of all 
``G_s`` values of the site. 


## Global radiation

Potential radiation for given time and latitude:
```@example doc
doy, hour = 160, 10.5
lat, long = 51.0, 11.5
potrad = potential_radiation(doy, hour, lat, long)
```

Calculation is based on sun's altitude, one of the horizontal coordinates of its position.
```@example doc
using Plots, StatsPlots, DataFrames, Dates, Pipe, Suppressor
hours = 0:24
lat,long = 51.0, 13.6 # Dresden Germany
#deg2second = 24*3600/360
doy = 160
datetimes = DateTime(2021) .+Day(doy-1) .+ Hour.(hours) #.- Second(round(long*deg2second))
res3 = @pipe calc_sun_position_hor.(datetimes, lat, long) |> DataFrame(_)
@df res3 scatter(datetimes, cols([:altitude,:azimuth]), legend = :topleft, # hide
  xlab="Date and Time", ylab = "rad", xrotation=6) # hide
```

The hour-angle at noon represents the difference to
local time. In the following example solar time is
about 55min ahead of local winter time.

```@example doc
summernoon = DateTime(2021) +Day(doy-1) + Hour(12) 
sunpos = calc_sun_position_hor(summernoon, lat, long) 
sunpos.hourangle * 24*60/(2*π) # convert angle to minutes
```

## Unit interconversions

The package further provides a number of useful unit interconversions, which are straightforward to use (please make sure that the input variable is in the right unit, e_g. rH has to be between 0 and 1 and not in percent):

```@example doc
# VPD to vapor pressure (e, kPa)
VPD_to_e(2.0, 25.0)
```
```@example doc
# vapor pressure to specific humidity (kg kg-1)
e_to_q(1.0, 100.0)
```
```@example doc
# relative humidity to VPD (kPa)
rH_to_VPD(0.6, 25.0)
```
```@example doc
# conductance from ms-1 to mol m-2 s-1
ms_to_mol(0.01, 25.0, 100.0) # mC, Tair, pressure
```
```@example doc
# umol CO2 m-2 s-1 to g C m-2 d-1
umolCO2_to_gC(20.0)
```

Many functions provide constant empirical parameters. Those can
be changed by overriding the default values with 
[`BigleafConstants`](@ref) 
and passing this Struct to the respective function.


