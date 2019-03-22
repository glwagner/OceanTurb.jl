var documenterSearchIndex = {"docs": [

{
    "location": "#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": ""
},

{
    "location": "#OceanTurb.jl-1",
    "page": "Home",
    "title": "OceanTurb.jl",
    "category": "section",
    "text": "OceanTurb.jl implements one-dimensional partial differential equations that model turbulent convection and diffusion.The primary purpose of the package is to explore models for convection and turbulent mixing in the ocean\'s surface boundary layer, where atmospheric forcing due to wind, waves, precipitation, evaporation, heating, cooling, and radiation drive turbulence and mediate the exchange of quantities like heat, momentum, and carbon between the atmosphere and ocean interior. Models that approximate the effects of atmospheric forcing on turbulence and turbulent mixing in the upper ocean are critical components in ocean circulation models and coupled climate models."
},

{
    "location": "#Installation-1",
    "page": "Home",
    "title": "Installation",
    "category": "section",
    "text": "Open julia, press ] at the Julian prompt to enter package manager mode, and typepkg> add https://github.com/glwagner/OceanTurb.jl.gitThis installs OceanTurb.jl from the master branch. The tests can then be run by typingpkg> test OceanTurbUse help mode by typing ? to find information about key functions:help?> iterate!\nsearch: iterate! iterate InteractiveUtils isinteractive Iterators\n\n  iterate!(model, Δt, nt=1)\n\n  Step model forward in time by one time-step with step-size Δt."
},

{
    "location": "#Components-1",
    "page": "Home",
    "title": "Components",
    "category": "section",
    "text": "Solvers for various turbulence models are implemented in submodules of OceanTurb.jl. For example, our simplest model solves the 1D diffusion equation and can be used by writingjulia> using OceanTurb.DiffusionIn addition to simple diffusion we haveK-Profile-Parameterization\n... others coming soon."
},

{
    "location": "#Authors-1",
    "page": "Home",
    "title": "Authors",
    "category": "section",
    "text": "Gregory L. Wagner."
},

{
    "location": "numerics/#",
    "page": "Numerical methods",
    "title": "Numerical methods",
    "category": "page",
    "text": ""
},

{
    "location": "numerics/#Numerical-methods-1",
    "page": "Numerical methods",
    "title": "Numerical methods",
    "category": "section",
    "text": "newcommandc \nnewcommandp \nnewcommanddpartial\n\nnewcommandr1mathrm1\n\nnewcommandeemathrme\n\nnewcommandbeqbeginequation\nnewcommandeeqendequation\n\nnewcommandbeqsbegingather\nnewcommandeeqsendgather"
},

{
    "location": "numerics/#Spatial-discretization-1",
    "page": "Numerical methods",
    "title": "Spatial discretization",
    "category": "section",
    "text": "OceanTurb.jl uses a one-dimensional finite-volume method to discretize momentum, temperature, salinity, and other variables.An ASCII-art respresentation of an example grid with nz=3 is ▲ z\n |              j=4   ===  Top   ▲              \n         i=3           *         | dzf (i=3)\n                j=3   ---        ▼\n         i=2           *             ▲            \n                j=2   ---            | dzc (j=2)\n         i=1           *             ▼  \n                j=1   ===  Bottomwhere the double lines indicate the top and bottom of the domain, the single lines indicate \"face\" boundaries, the i\'s index cell centers (nodes) and j\'s index the z-location of cell interfaces (faces). Horizontal momentum and tracer variables are located at cell centers, while fluxes of these quantities (and vertical-velocity-like variables when present) are located at cell faces."
},

{
    "location": "numerics/#Derivatives-and-diffusive-fluxes-1",
    "page": "Numerical methods",
    "title": "Derivatives and diffusive fluxes",
    "category": "section",
    "text": "The derivative of a quantity phi at face i isbeq\nleft( d_z phi right )_i = fracphi_i - phi_i-1Delta c_i c\neeqwhere phi_i denotes the value of phi at cell i, and Delta c_i = z_c i - z_c i-1 is the vertical separation between node i and cell point i-1.With diffusivity given on cell interfaces, the diffusive flux at face i is just K_i left ( d_z phi right )_i. The divergence of the diffusive flux at node i is thereforenewcommandKdz1K_1 left ( d_z phi right )_1 \nbeginalign\nleft ( d_z K d_z phi right )_i = frac Kdzi+1 - Kdzi Delta f_i c \n= frac\n          tfracK_i+1Delta c_i+1 phi_i+1\n        - left ( tfracK_i+1Delta c_i+1 + tfracK_iDelta c_i right ) phi_i\n         + tfracK_iDelta c_i phi_i-1Delta f_i\nendalign"
},

{
    "location": "numerics/#Time-stepping-1",
    "page": "Numerical methods",
    "title": "Time-stepping",
    "category": "section",
    "text": "To integrate ocean surface boundary layer models forward in time, we implement explicit time-stepping schemes for partial differential equations of the formbeq\nd_t Phi = R(Phi) c\nlabelexplicitform\neeqwhere Phi is a variable like velocity, temperature, or salinity and R is an arbitrary function representing any number of processes, including turbulent diffusion and internal forcing.In the future, we will implement implicit-explicit schemes for partial differential equations of the formbeq labelimplicitdiffusion\nd_t Phi - d_z K d_z phi = R(Phi) c\neeqwhere K is a diffusivity. These implicit-explicit schemes will treat diffusive terms on the left of \\eqref{implicitdiffusion} implicitly."
},

{
    "location": "numerics/#Forward-Euler-method-1",
    "page": "Numerical methods",
    "title": "Forward-Euler method",
    "category": "section",
    "text": "An explicit forward Euler time integration scheme uses the temporal discretization \\eqref{explicitform}:beq\nPhi^n+1 = Phi^n + Delta t R left ( Phi^n right )\neeqwhere the superscripts n and n+1 denote the solution at time-step n and n+1, respectively."
},

{
    "location": "models/basics/#",
    "page": "Primer on boundary layer modeling",
    "title": "Primer on boundary layer modeling",
    "category": "page",
    "text": ""
},

{
    "location": "models/basics/#Primer-on-boundary-layer-modeling-1",
    "page": "Primer on boundary layer modeling",
    "title": "Primer on boundary layer modeling",
    "category": "section",
    "text": "newcommandc \nnewcommandp \nnewcommanddpartial\n\nnewcommandr1mathrm1\nnewcommandb1boldsymbol1\n\nnewcommandeemathrme\n\nnewcommandbeqbeginequation\nnewcommandeeqendequation\n\nnewcommandbeqsbegingather\nnewcommandeeqsendgatherModels for the ocean surface boundary layer are partial differential equations that approximate the effects of atmospheric forcing on the turbulent vertical flux and evolution of large-scale temperature, salinity, and momentum fields.Internal and surface fluxes of heat are due toabsorption of incoming shortwave solar radiation;\ncooling by outgoing longwave radiation;\nlatent and sensible heat exchange with the atmosphere.Surface fluxes of salinity occur due to evaporation and precipitation, while momentum fluxes are associated with atmospheric winds.Vertical turbulent fluxes are typically associated withgravitational instability and convection, and\nmechanical turbulent mixing associated with currents and wind forcing.OceanTurb.jl uses an implementation of atmospheric and radiative forcings that is shared across all models. The models therefore differ in the way they parameterize convective and wind-driven mechanical mixing."
},

{
    "location": "models/basics/#Coordinate-system-1",
    "page": "Primer on boundary layer modeling",
    "title": "Coordinate system",
    "category": "section",
    "text": "We use a Cartesian coordinate system in which gravity points downwards, toward the ground or bottom of the ocean. The vertical coordinate z thus increases upwards. We locate the surface at z=0. This means that if the boundary layer has depth h, the bottom of the boundary layer is located at z=-h."
},

{
    "location": "models/basics/#Governing-equations-1",
    "page": "Primer on boundary layer modeling",
    "title": "Governing equations",
    "category": "section",
    "text": "The one-dimensional, horizontally-averaged boundary-layer equations for horizontal momentum U and V, salinity S, and temperature T arebeqs\nU_t =   f V - d_z overlinew u      c labelxmomentum \nV_t = - f U - d_z overlinew v      c \nT_t =       - d_z overlinew theta + I_theta c labeltemperature \nS_t =       - d_z overlinew s      c labelsalinity \neeqswhere subscripts t and z denote derivatives with respect to time and the vertical coordinate z and f is the Coriolis parameter. The lowercase variables u, v, s, and theta refer to the three-dimensional perturbations from horizontal velocity, salinity, and temperature, respectively. In \\eqref{xmomentum}–\\eqref{temperature}, internal forcing of temperature due to solar radiation is denoted I_Phi."
},

{
    "location": "models/basics/#Surface-fluxes-1",
    "page": "Primer on boundary layer modeling",
    "title": "Surface fluxes",
    "category": "section",
    "text": "Turbulence in the ocean surface boundary layer is driven by fluxes from the atmosphere above. A surface flux of some variable phi is denoted F_phi. Surface fluxes includeMomentum fluxes due to wind, denoted F_u bx + F_v by = -rho_0 btau for wind stress btau;\nTemperature flux F_theta = - Q  rho_0 c_P associated with \'heating\' Q;\nSalinity flux F_s = (E-P)S associated evaporation E and precipitation P.We use the traditional convention ordinary to physics, but not always ordinary to oceanography, in which a  positive flux corresponds to the movement of a quantity in the positive z-direction. This means, for example, that a positive vertical velocity w gives rise to a positive advective flux w phi. This convention also implies that a positive temperature flux at the ocean surface –- corresponding to heat fluxing upwards, out of the ocean, into the atmosphere –- implies a cooling of the ocean surface boundary layer."
},

{
    "location": "models/kpp/#",
    "page": "The K-Profile-Parameterization (KPP)",
    "title": "The K-Profile-Parameterization (KPP)",
    "category": "page",
    "text": ""
},

{
    "location": "models/kpp/#The-K-Profile-Parameterization-(KPP)-1",
    "page": "The K-Profile-Parameterization (KPP)",
    "title": "The K-Profile-Parameterization (KPP)",
    "category": "section",
    "text": "newcommandc          \nnewcommandp          \nnewcommandd         partial\nnewcommandr1      mathrm1\nnewcommandb1      boldsymbol1\nnewcommandee        mathrme\nnewcommanddi         mathrmd\nnewcommandep        epsilon\n\nnewcommandbeq       beginequation\nnewcommandeeq       endequation\nnewcommandbeqs      begingather\nnewcommandeeqs      endgather\n\n Non-dimensional numbers\nnewcommandRi        mathrmRi\nnewcommandK         mathrmKE        \n\nnewcommandbtau      btau  wind stress vector\n\n Model functions and constants\nrenewcommandF2      Upsilon^1_2\nrenewcommandC2      C^1_2\n\nnewcommanduwind     varpi_tau\nnewcommandubuoy     varpi_bThe K-Profile-Parameterization, or \"KPP\", was proposed by Large et al (1994) as a model for convection- and wind-driven mixing in the upper ocean. In KPP, vertical turbulent fluxes of a quantity phi are parameterized asbeq\noverlinew phi = - K_phi d_z Phi + N_phi c\neeqwhere Phi is the resolved or horizontally-averaged quantity, K_phi is a turbulent diffusivity, and N_phi is a \'non-local\' flux.The non-local flux and turbulent diffusivity are defined to vanish at the surface, and at the bottom of the \'mixing layer\' h, which roughly corresponds to the depth at which turbulent fluxes and turbulent kinetic energy decay to zero."
},

{
    "location": "models/kpp/#Description-of-the-model-1",
    "page": "The K-Profile-Parameterization (KPP)",
    "title": "Description of the model",
    "category": "section",
    "text": "Below, we denote \'model parameters\' as Cmathrmlabelmathrmvar, where \'label\' describes the parameter and \'var\' is a variable like U V T or S.Buoyancy is then defined in terms of T and S asbeginalign\nB  equiv - fracg rhorho_0 \n       = g left  alpha left ( T - T_0 right ) - beta left ( S - S_0 right ) right  c\nendalignwhere g = 981  mathrmm  s^-2 alpha = 2 times 10^-4  mathrmK^-1, and beta = 8 times 10^-5 mathrmK^-1, are gravitational acceleration, the thermal expansion coefficient, and the haline contraction coefficient, respectively. Buoyancy forcing is equivalentlybeq\nF_b = g left ( alpha F_theta - beta F_s right ) p\neeqThe turbulent velocity scales associated with buoyancy and wind forcing arebeq\nomega_b equiv  h F_b ^13 qquad textand qquad omega_tau equiv  btau  rho_0 ^12 p\neeqwhere h is the depth of the \'mixing layer\', or the depth to which turbulent mixing and turbulent fluxes penetrate, btau is wind stress, and rho_0 = 1028  mathrmkg  m^-3 is a reference density.We also define the ratiosbeq\nr_b equiv left ( fracomega_bomega_tau right )^3 qquad textand\nqquad r_tau equiv left ( fracomega_tauomega_b right )^3 = frac1r_b p\neeq"
},

{
    "location": "models/kpp/#The-boundary-layer-depth,-h-1",
    "page": "The K-Profile-Parameterization (KPP)",
    "title": "The boundary layer depth, h",
    "category": "section",
    "text": "The boundary layer depth h is defined implicitly via the bulk Richardson number criterionbeq labelbulk_ri\nCRi = frach Delta B(-h) Delta bU(-h)^2 + FK(-h) p\neeqwhere the critical Ri is CRi = 03. The operator Delta is definedbeq\nDelta Phi(z) = -frac1Cep z int_Cep z^0 Phi(z) di z - Phi(z) c\neeqwhere Cep = 01 is the surface layer fraction. The function FK(z) isbeq  labelunresolved_ke\nFK(z) = CK (-z)^43 sqrt max left  0 B_z(z) right   max left  0 F_b^13 right  p\neeqThe unresolved kinetic energy constant is CK = 432.where g is gravitational acceleration and rho_0 T_0 S_0 are reference densities, temperatures, and salinities. rho is the density deviation from the reference."
},

{
    "location": "models/kpp/#Non-local-flux-1",
    "page": "The K-Profile-Parameterization (KPP)",
    "title": "Non-local flux",
    "category": "section",
    "text": "The non-local flux is defined only for T and S, and isbeq\nN_Phi = CN F_Phi d (1 - d)^2 c\neeqwhere d = -zh is a non-dimensional depth coordinate and CN = 633."
},

{
    "location": "models/kpp/#Turbulent-Diffusivity-1",
    "page": "The K-Profile-Parameterization (KPP)",
    "title": "Turbulent Diffusivity",
    "category": "section",
    "text": "The KPP diffusivity is definedbeq\nK_phi = h FwPhi d ( 1 - d )^2 c\neeqwhere FwPhi is the turbulent velocity scale. In wind-driven turbulence under stable buoyancy forcing such that F_b  0, the turbulent velocity scale isbeq\nFwPhi = frac Ckappa uwind1 + Cmathrmstab r_b d p\neeqwherea Ckappa = 04 is the Von Karman constant and Cmathrmstab = 20.In wind-driven turbulence but under destabilizing buoyancy forcing, when min left  Cep d right   Cdphi r_tau, the turbulent velocity scale isbeq\nFwPhi = Ckappa omega_tau left ( 1 + Cmathrmunst r_b min left  Cep d right  right )^n_Phi c\neeqwhere n_U = 14 for velocities, n_T = 12 for tracers, Cmathrmunst = 64, the transition parameter for velocities is CdU = 05, and the transition parameter for tracers is CdT = 25.In convection-driven turbulence affected by wind mixing, when min left  Cep d right  = Cdphi r_tau, the turbulent velocity scale isbeq\nFwPhi = CbPhi omega_b left ( min left  Cep d right  + CtauPhi r_tau right )^13 c\neeqwhere CbU = CbV = 0215, CbT = CbS = 253, CtauU = CtauV = 00806, and CtauT = CtauS = 185."
},

{
    "location": "models/kpp/#Selected-tests-1",
    "page": "The K-Profile-Parameterization (KPP)",
    "title": "Selected tests",
    "category": "section",
    "text": "See /test/runtests.jl for more tests."
},

{
    "location": "models/kpp/#Linear-temperature-profile-and-no-velocity-field-1",
    "page": "The K-Profile-Parameterization (KPP)",
    "title": "Linear temperature profile and no velocity field",
    "category": "section",
    "text": "Some simple tests can be defined when the model state is U=V=S=0 and T = gamma z. In this case the buoyancy becomes B = g alpha gamma z and the buoyancy gradient is B_z = g alpha gamma. If we further take Cep to 0 and F_b  0, and note that the value of T in the top grid cell is T_N = -gamma Delta z  2, where N is the number of grid points, we find thatbeq\nDelta T(-h) = gamma h - T_N c\neeqandbeq\nDelta B(-h) = g alpha gamma h - B_N p\neeqThe unresolved kinetic energy function isbeq\nFK(-h) = CK h^43 sqrtg alpha gamma F_b^13 p\neeqThe bulk Richardson number criterion then becomesbeginalign\nCRi = frach Delta B(-h)FK(-h) c \n          = fracg alpha gamma h^2 - h B_NCK h^43 sqrtg alpha gamma F_b^13 c \n          = fracg alpha gamma h - B_NCK sqrtg alpha gamma left ( h F_b right )^13 c \nendalignModifying the temperature profile so that T_N = B_N = 0 allows us to analytically calculate the mixed layer depth:beq\nh = left ( CRi CK right )^32 sqrtF_b left ( g alpha gamma right )^-34 p\nlabelanalyticaldepth\neeq"
},

{
    "location": "models/kpp/#Linear-temperature-profile-and-linear-shear-1",
    "page": "The K-Profile-Parameterization (KPP)",
    "title": "Linear temperature profile and linear shear",
    "category": "section",
    "text": "Consider a piecewise-constant temperature profilebeq\nT(z) =  left  beginmatrix\n  T_0  quad textfor  z  -h c \n  -T_0  quad textfor  z  -h c\n  endmatrix right \neeqand a velocity profilebeq\nU(z) =  left  beginmatrix\n  U_0  quad textfor  z  -h c \n  -U_0  quad textfor  z  -h c\n  endmatrix right \neeqso that T(z=-h) = U(z=-h) = 0. We then have Delta T(-h) = T_0 and Delta U^2(-h) = U_0^2, so that with g=alpha=1,beq\nCRi = frach T_0U_0^2 p\nlabelsheardepth\neeqSetting h = 9, CRi=1, T_0 = 1, and U_0=3 yields a consistent solution."
},

{
    "location": "models/kpp/#Limiting-cases-for-turbulent-velocity-scales-1",
    "page": "The K-Profile-Parameterization (KPP)",
    "title": "Limiting cases for turbulent velocity scales",
    "category": "section",
    "text": "Under zero momentum forcing, the turbulent vertical velocity scale isbeq\nFwPhi = CbPhi left ( Cep right )^13  h F_b ^13 p\nlabelbuoyscaletest\neeqwe write the test in \\eqref{buoyscaletest} using the depth in \\eqref{analyticaldepth}Under zero buoyancy forcing, the turbulent velocity scale isbeq\nFwPhi = Ckappa omega_tau p\nlabelwindscaletest\neeq"
},

{
    "location": "models/kpp/#Table-of-model-parameters-1",
    "page": "The K-Profile-Parameterization (KPP)",
    "title": "Table of model parameters",
    "category": "section",
    "text": "The model parameters in KPP areParameter Value Reference\nCep 0.1 pretty much everywhere"
},

{
    "location": "models/kpp/#References-1",
    "page": "The K-Profile-Parameterization (KPP)",
    "title": "References",
    "category": "section",
    "text": "Large et al (1994)\nCVMix documentation\nVan Roekel et al (2018)"
},

{
    "location": "models/pacanowskiphilander/#",
    "page": "Pacanowski-Philander",
    "title": "Pacanowski-Philander",
    "category": "page",
    "text": ""
},

{
    "location": "models/pacanowskiphilander/#Pacanowski-Philander-1",
    "page": "Pacanowski-Philander",
    "title": "Pacanowski-Philander",
    "category": "section",
    "text": "newcommandc \nnewcommandp \nnewcommanddpartial\n\nnewcommandr1mathrm1\n\nnewcommandeemathrme\n\nnewcommandbeqbeginequation\nnewcommandeeqendequation\n\nnewcommandbeqsbegingather\nnewcommandeeqsendgatherPacanowski and Philander (1981) propose a simple one-dimensional model for equatorial boundary layers dominated by mechanical turbulent mixing. In their model, horizontal velocity, temperature, and salinity are governed bybeqs\nU_t =   f V - d_z overlinew u - F^u c \nV_t = - f U - d_z overlinew v - F^v c \nS_t =       - d_z overlinew s - F^S c \nT_t =       - d_z overlinew theta - F^T c\neeqswhere uppercase variables are resolved, mean quantities, and lowercase variables are unresolved perturbations. A key variable in the Pacaonwski and Philander formulation is the local Richardson number, defined byRi = - fracg rho_zrho_0 left ( U_z^2 + V_z^2 right ) cwhere g is the gravitational constant, rho_0 is a reference density, and rho is density. For simplicity, we use a linear equation of state between temperature, salinity, and density,rho(z t) = rho_0 left \n  alpha left ( T-T_0 right ) + beta left ( S-S_0 right ) right  cwhere typical values for the temperature and salinity expansion coefficients are alpha = 2 times 10^-4 rK m^3  kg and beta = 1."
},

{
    "location": "models/pacanowskiphilander/#Eddy-diffusivities-for-momentum,-temperature,-and-salinity-1",
    "page": "Pacanowski-Philander",
    "title": "Eddy diffusivities for momentum, temperature, and salinity",
    "category": "section",
    "text": "The core of the PP81 model is to parameterize vertical turbulent fluxes with an eddy diffusivity/eddy viscosity, such thatoverlinew phi = kappa_Phi d_z Phi cwhere kappa_Phi is the eddy diffusivity or viscosity of the quantity Phi.The eddy viscosity for both U and V iskappa_U = nu_0 + fracnu_1left ( 1 + c Ri right )^n cwhile the eddy diffusivity for T and S arekappa_T = kappa_0 + frac kappa_1  left ( 1 + c Ri right )^n+1 pThis parameterization implies an Ri-dependent turbulent Prandtl number ofPr = frackappa_Ukappa_T approx left ( 1 + c Ri right ) cwhen Ri is small. Typical values for the parameters (see CV12 arenu_0 = 10^-4 rcms,\nnu_1 = 10^-2 rcms,\nkappa_0 = 10^-5 rcms,\nkappa_1 = 10^-2 rcms,\nc = 5, and\nn=2."
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
    "text": "Modules = [OceanTurb]\nPublic = false\nOrder = [:type]"
},

{
    "location": "man/functions/#",
    "page": "Functions",
    "title": "Functions",
    "category": "page",
    "text": ""
},

{
    "location": "man/functions/#Base.Libc.time-Tuple{AbstractModel}",
    "page": "Functions",
    "title": "Base.Libc.time",
    "category": "method",
    "text": "Get the current simulation time of the model.\n\n\n\n\n\n"
},

{
    "location": "man/functions/#OceanTurb.FluxBoundaryCondition-Tuple{Any}",
    "page": "Functions",
    "title": "OceanTurb.FluxBoundaryCondition",
    "category": "method",
    "text": "FluxBoundaryCondition(boundary, flux)\n\nConstuct a flux boundary condition that specifies the flux of some field on a boundary. If flux is a function, its arguments must be synced with the expection of Model.\n\n\n\n\n\n"
},

{
    "location": "man/functions/#OceanTurb.ValueBoundaryCondition-Tuple{Any}",
    "page": "Functions",
    "title": "OceanTurb.ValueBoundaryCondition",
    "category": "method",
    "text": "ValueBoundaryCondition(boundary, flux)\n\nConstuct a flux boundary condition that specifies the flux of some field on a boundary. If flux is a function, its arguments must be synced with the expection of Model.\n\n\n\n\n\n"
},

{
    "location": "man/functions/#OceanTurb.ZeroFluxBoundaryConditions-Tuple{}",
    "page": "Functions",
    "title": "OceanTurb.ZeroFluxBoundaryConditions",
    "category": "method",
    "text": "Returns FieldBoundaryConditions with zero flux at top and bottom.\n\n\n\n\n\n"
},

{
    "location": "man/functions/#OceanTurb.iter-Tuple{AbstractModel}",
    "page": "Functions",
    "title": "OceanTurb.iter",
    "category": "method",
    "text": "Get the current iteration of the model.\n\n\n\n\n\n"
},

{
    "location": "man/functions/#OceanTurb.iterate!-Tuple{Any,Any,Any}",
    "page": "Functions",
    "title": "OceanTurb.iterate!",
    "category": "method",
    "text": "iterate!(model, Δt, nt=1)\n\nStep model forward in time by one time-step with step-size Δt.\n\n\n\n\n\n"
},

{
    "location": "man/functions/#OceanTurb.set!-Tuple{AbstractSolution}",
    "page": "Functions",
    "title": "OceanTurb.set!",
    "category": "method",
    "text": "set!(solution, kwargs...)\n\nSet the fields of a solution. For example, use\n\nT0 = rand(4) S0(z) = exp(-z^2/10) set!(solution, T=T0, S=S0)\n\nTo set solution.T and solution.S to T0 and S0.\n\n\n\n\n\n"
},

{
    "location": "man/functions/#OceanTurb.set_bcs!-Tuple{Any,Any,Tuple}",
    "page": "Functions",
    "title": "OceanTurb.set_bcs!",
    "category": "method",
    "text": "set_bcs!(model, fld, bcs::Tuple)\n\nSet boundary conditions on fld as  bcs = (bottom_bc, top_bc)\n\nset_bcs!(model; kwargs...)\n\nSet boundary conditions, where the keywords are fields of model.solution and their arguments are tuples of (bottom, top) boundary condition.\n\n\n\n\n\n"
},

{
    "location": "man/functions/#OceanTurb.set_bottom_bc!-Tuple{Any,Any,Any}",
    "page": "Functions",
    "title": "OceanTurb.set_bottom_bc!",
    "category": "method",
    "text": "Set the bottom boundary condition for fld in model.\n\n\n\n\n\n"
},

{
    "location": "man/functions/#OceanTurb.set_flux_bcs!-Tuple{Any}",
    "page": "Functions",
    "title": "OceanTurb.set_flux_bcs!",
    "category": "method",
    "text": "set_flux_bcs!(model, fld, bcs::Tuple)\n\nSet flux boundary conditions on fld as  bcs = (bottom_bc, top_bc)\n\nset_flux_bcs!(model; kwargs...)\n\nSet flux boundary conditions for model.\n\n\n\n\n\n"
},

{
    "location": "man/functions/#OceanTurb.set_top_bc!-Tuple{Any,Any,Any}",
    "page": "Functions",
    "title": "OceanTurb.set_top_bc!",
    "category": "method",
    "text": "Set the top boundary condition for fld in model.\n\n\n\n\n\n"
},

{
    "location": "man/functions/#OceanTurb.Δc-Tuple{UniformGrid,Any}",
    "page": "Functions",
    "title": "OceanTurb.Δc",
    "category": "method",
    "text": "Return the cell spacing at index i.\n\n\n\n\n\n"
},

{
    "location": "man/functions/#OceanTurb.Δf-Tuple{UniformGrid,Any}",
    "page": "Functions",
    "title": "OceanTurb.Δf",
    "category": "method",
    "text": "Return the face spacing at index i.\n\n\n\n\n\n"
},

{
    "location": "man/functions/#OceanTurb.∂z!-Tuple{Field{Cell,A,G} where G where A,Field{Face,A,G} where G where A}",
    "page": "Functions",
    "title": "OceanTurb.∂z!",
    "category": "method",
    "text": "Calculate c = ∂f/∂z in the grid interior.\n\n\n\n\n\n"
},

{
    "location": "man/functions/#OceanTurb.∂z!-Tuple{Field{Face,A,G} where G where A,Field{Cell,A,G} where G where A}",
    "page": "Functions",
    "title": "OceanTurb.∂z!",
    "category": "method",
    "text": "Calculate f = ∂c/∂z in the grid interior.\n\n\n\n\n\n"
},

{
    "location": "man/functions/#OceanTurb.∂z-Tuple{Any,Any}",
    "page": "Functions",
    "title": "OceanTurb.∂z",
    "category": "method",
    "text": "∂z(a, i)\n\nReturn the discrete derivative of a at grid point i.\n\nThe derivative of a Field{Cell} is computed at face points, and the derviative of a Field{Face} is computed at cell points.\n\n\n\n\n\n"
},

{
    "location": "man/functions/#OceanTurb.∂z-Tuple{Any}",
    "page": "Functions",
    "title": "OceanTurb.∂z",
    "category": "method",
    "text": "Return the CellField ∂f/∂z, where f is a FaceField.\n\n\n\n\n\n"
},

{
    "location": "man/functions/#OceanTurb.∂z-Tuple{Field{Cell,A,G} where G where A,Any}",
    "page": "Functions",
    "title": "OceanTurb.∂z",
    "category": "method",
    "text": "Return ∂c/∂z at face index i.\n\n\n\n\n\n"
},

{
    "location": "man/functions/#OceanTurb.∂z-Tuple{Field{Cell,A,G} where G where A}",
    "page": "Functions",
    "title": "OceanTurb.∂z",
    "category": "method",
    "text": "Return the FaceField ∂c/∂z, where c is a CellField.\n\n\n\n\n\n"
},

{
    "location": "man/functions/#OceanTurb.∂z-Tuple{Field{Face,A,G} where G where A,Any}",
    "page": "Functions",
    "title": "OceanTurb.∂z",
    "category": "method",
    "text": "Return ∂c/∂z at face index i.\n\n\n\n\n\n"
},

{
    "location": "man/functions/#Functions-1",
    "page": "Functions",
    "title": "Functions",
    "category": "section",
    "text": "Modules = [OceanTurb]\nPrivate = false\nOrder = [:function]"
},

]}
