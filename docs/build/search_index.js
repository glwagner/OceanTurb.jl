var documenterSearchIndex = {"docs": [

{
    "location": "#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": "newcommandc \n\nnewcommandr1mathrm1\n\nnewcommandeemathrme\n\nnewcommandbeqbeginequation\nnewcommandeeqendequation\n\nnewcommandbeqsbegingather\nnewcommandeeqsendgather"
},

{
    "location": "#OceanBoundaryLayerModels.jl-1",
    "page": "Home",
    "title": "OceanBoundaryLayerModels.jl",
    "category": "section",
    "text": "OceanMixedLayers.jl implements models for the ocean\'s boundary layer."
},

{
    "location": "#Physics-of-the-oceanic-boundary-layer-1",
    "page": "Home",
    "title": "Physics of the oceanic boundary layer",
    "category": "section",
    "text": "The dynamics of the ocean\'s boundary layer are dominated by forcing from the atmosphere. Atmospheric forcing includes momentum forcing by winds, salinity forcing by evaporation and precipitation, and heat forcing by latent heat fluxes, sensible heat fluxes, and incoming and outgoing radiation."
},

{
    "location": "#Governing-equations-for-1D-boundary-layers-1",
    "page": "Home",
    "title": "Governing equations for 1D boundary layers",
    "category": "section",
    "text": "We write the one-dimensional boundary-layer balances for horizontal momentum,  salinity, and temperature asbeqs\nu_t - f v = -G^u_z - F^u_z c \nv_t + f u = -G^v_z - F^v_z c \n      S_t = -G^S_z - F^S_z c \n      T_t = -G^T_z - F^T_z c\neeqswhere G^phi = overlinew phi denotes the turbulent vertical flux of  a quantity phi, while F^phi denotes vertical fluxes due to  forcing."
},

{
    "location": "#Temperature-forcing-1",
    "page": "Home",
    "title": "Temperature forcing",
    "category": "section",
    "text": "We write the temperature forcing F^T asF^T = F^rlat + F^rsens + F^rlongwave \n        + F^rshortwave cin terms of the four contributions from latent heating, sensible heating,  outgoing longwave radiation, and incoming shortwave radiation. The first three contributions are implemented as boundary conditions. Shortwave  radiation, on the other hand, heats the interior of the boundary layer. We  parameterize the effect of interior heating by shortwave radiation by dividing the specrtum into infrared (IR) and ultraviolet (UV) components and introducing an absorption functionbeq\nF^rshortwave = F^rshortwave_0\n    left ( alpha_rIR exp left  zd_rIR right  \n          + alpha_rUV exp left  zd_rUV right  right ) c\nlabelshortwaverad\neeqwhere F^rshortwave_0 is the incoming shortwave radiation at the  surface, which is provided as an input to the boundary layer model. In \\eqref{shortwaverad},  d_rIR and d_rUV are the  penetration scales of infrared and ultraviolet radiation, and alpha_rIR and alpha_rUV are the fractions of total shortwave radiation with infrared and ultraviolet wavelengths, respectively, such that  alpha_rIR + alpha_rUV = 1."
},

{
    "location": "man/types/#",
    "page": "Private types",
    "title": "Private types",
    "category": "page",
    "text": ""
},

{
    "location": "man/types/#Private-types-1",
    "page": "Private types",
    "title": "Private types",
    "category": "section",
    "text": ""
},

{
    "location": "man/types/#Private-types-in-module-OceanMixedLayerModels:-1",
    "page": "Private types",
    "title": "Private types in module OceanMixedLayerModels:",
    "category": "section",
    "text": "Modules = [OceanMixedLayerModels]\nPublic = false\nOrder = [:type]"
},

{
    "location": "man/functions/#",
    "page": "Functions",
    "title": "Functions",
    "category": "page",
    "text": ""
},

{
    "location": "man/functions/#Functions-1",
    "page": "Functions",
    "title": "Functions",
    "category": "section",
    "text": ""
},

{
    "location": "man/functions/#OceanMixedLayerModels.loadforcing-Tuple{Any}",
    "page": "Functions",
    "title": "OceanMixedLayerModels.loadforcing",
    "category": "method",
    "text": "loadforcing(filepath)\n\nInitalize a Forcing from the JLD2 file at filepath. The names of the file\'s  fields must correspond to the names of the Forcing fields.\n\n\n\n\n\n"
},

{
    "location": "man/functions/#Functions-exported-from-OceanMixedLayerModels:-1",
    "page": "Functions",
    "title": "Functions exported from OceanMixedLayerModels:",
    "category": "section",
    "text": "Modules = [OceanMixedLayerModels]\nPrivate = false\nOrder = [:function]"
},

]}
