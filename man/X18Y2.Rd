\name{X18Y2}
\alias{X18Y2}
\docType{data}
\title{
Aircraft InfraRed Signature (IRS)
}
\description{
  Aircraft flight measures. The eight first ones
are related to flight conditions and IR optical
  properties of aircraft surfaces.
  The  eight following ones are related to atmospheric conditions.
  The last two variables are the response values.
Three qualitative measures have been transformed into 
numeric variables: MODEL and CLOUDS are coded as 0/1
  and IHAZE is split into three 0/1 variables.
}
\usage{data("X18Y2")}
\format{
  A data frame with 180 observations on the following 20 variables.
  All are numeric.
 \tabular{ll}{
    ALTI \tab  Altitude\cr
    MACH \tab  Mach number\cr
    POWERS \tab  Engine power setting\cr
    EAI \tab  Emissivity of air intake\cr
    CAP \tab  Cap\cr
    YAW \tab  Yaw\cr
    ROLL \tab  Roll\cr
    PITCH \tab  Pitch angles\cr
    VIS \tab  Visibility\cr
    RH \tab  Relative humidity\cr
    TA \tab  Ground air temperature\cr
    MODEL \tab  Atmospheric model\cr
    IHAZE1 \tab  Aerosol model urban\cr
    IHAZE2 \tab  Aerosol model maritime\cr
    IHAZE3 \tab  Aerosol mode rural\cr
    CLOUDS \tab  Cloud presence\cr
    HBASE \tab Base altitude of
the cloud layer\cr
    HOUR \tab  Hour to compute solar position\cr
    IRStot \tab  Typical infrared sensor spectral ranges\cr
    IRMS \tab  Multispectral intensity\cr
  }
}
\source{
Lefebvre, S. ONERA, Palaiseau, France.
}
\references{
 Gauchi, J.-P. 2015.
      A practical method of global sensitivity analysis
under constraints.
Rapport technique 2015-1.
INRA, UR1404, F-78350 Jouy-en-Josas, France.
}
\keyword{datasets}
