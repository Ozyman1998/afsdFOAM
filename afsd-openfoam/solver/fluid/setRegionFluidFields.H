    // Referencias a los campos de la región fluida i
    volScalarField& alpha1 = alpha1Fluid[i];
    volScalarField& alpha2 = alpha2Fluid[i];
    volVectorField& U = UFluid[i];
    volScalarField& p = pFluid[i];
    volScalarField& T = TFluid[i];

    surfaceScalarField& phi = phiFluid[i];
    surfaceScalarField& alphaPhi = alphaPhiFluid[i];
    surfaceScalarField& rhoPhi = rhoPhiFluid[i];
    surfaceScalarField& rhoCpPhi = rhoCpPhiFluid[i];

    volScalarField& rho = rhoFluid[i];
    volScalarField& rhoCp = rhoCpFluid[i];
    volScalarField& mu = muFluid[i];
    volScalarField& muMetal = muMetalFluid[i];
    volScalarField& kappa = kappaFluid[i];
    volScalarField& epsDot = epsDotFluid[i];
    volScalarField& Qvisc = QviscFluid[i];

    fv::options& fvOptions = fluidFvOptions[i];

    // Herramienta (tesis 2.2.1.1): constant/<region>/toolProperties
    const dictionary& toolDict = toolPropsFluid[i];

    const label pRefCell = pRefCellFluid[i];
    const scalar pRefValue = pRefValueFluid[i];

    // Propiedades de las fases (tesis cap. 2, ec. 2 y 11)
    const dictionary& phaseProps = phasePropsFluid[i];
    const dictionary& metalDict = phaseProps.subDict("metal");
    const dictionary& airDict = phaseProps.subDict("air");

    const dimensionedScalar rho1("rho", dimDensity, metalDict);
    const dimensionedScalar rho2("rho", dimDensity, airDict);
    const dimensionedScalar Cp1("Cp", dimSpecificHeatCapacity, metalDict);
    const dimensionedScalar Cp2("Cp", dimSpecificHeatCapacity, airDict);
    const dimensionedScalar mu1
    (
        "mu",
        dimDynamicViscosity,
        metalDict.getOrDefault<scalar>("mu", 0)
    );
    const dimensionedScalar mu2("mu", dimDynamicViscosity, airDict);
    const scalarList kappaCoeffs1(metalDict.get<scalarList>("kappaCoeffs"));
    const scalarList kappaCoeffs2(airDict.get<scalarList>("kappaCoeffs"));
    const scalar betaTQ(phaseProps.get<scalar>("betaTQ"));

    // Rango de validez del ajuste k(T) [°C]
    const scalar kappaTmin(phaseProps.getOrDefault<scalar>("kappaTmin", 0));
    const scalar kappaTmax(phaseProps.getOrDefault<scalar>("kappaTmax", 600));

    // Viscosidad del metal: "constant" (mu) o "SellarsTegart" (ec. 6-9)
    const word metalViscosityModel
    (
        metalDict.getOrDefault<word>("viscosityModel", "constant")
    );
    const bool isSellarsTegart = (metalViscosityModel == "SellarsTegart");
    const dictionary stDict
    (
        metalDict.subOrEmptyDict("SellarsTegartCoeffs")
    );

    // Tabla 2.1
    const scalar stSigmaN(isSellarsTegart ? stDict.get<scalar>("sigmaN") : 0);
    const scalar stA(isSellarsTegart ? stDict.get<scalar>("A") : 1);
    const scalar stN(isSellarsTegart ? stDict.get<scalar>("n") : 1);
    const scalar stQ(isSellarsTegart ? stDict.get<scalar>("Q") : 0);
    const scalar stR(isSellarsTegart ? stDict.get<scalar>("R") : 1);

    // Límites numéricos (no están en la tesis)
    const scalar stEpsDotMin(isSellarsTegart ? stDict.get<scalar>("epsDotMin") : 1);
    const scalar stMuMax(isSellarsTegart ? stDict.get<scalar>("muMax") : GREAT);
    const scalar stTmin(isSellarsTegart ? stDict.get<scalar>("Tmin") : 0);
    const scalar stTmax(isSellarsTegart ? stDict.get<scalar>("Tmax") : GREAT);
