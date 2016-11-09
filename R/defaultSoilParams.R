defaultSoilParams<-function() {
  return(list(
    SoilDepth = 1000,
    RockLayerDepth = 4000,
    TS_clay = 25,
    TS_silt = 50,
    TS_sand = 25,
    TS_macro = 0.1,
    TS_rfc = 20,
    SS_clay = 25,
    SS_silt = 50,
    SS_sand = 25, 
    SS_macro = 0.1,
    SS_rfc = 40,
    RL_clay = 25,
    RL_silt = 50,
    RL_sand = 25, 
    RL_macro = 0.1,
    RL_rfc = 85,
    Gsoil= 0.5,
    Ksoil=0.05));
}