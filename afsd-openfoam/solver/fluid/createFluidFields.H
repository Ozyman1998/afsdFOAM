// Región build: mezcla inmiscible metal (fase 1) + aire (fase 2), VOF.
// Propiedades en constant/<region>/phaseProperties.

PtrList<IOdictionary> phasePropsFluid(fluidRegions.size());
PtrList<IOdictionary> toolPropsFluid(fluidRegions.size());

PtrList<volScalarField> alpha1Fluid(fluidRegions.size());
PtrList<volScalarField> alpha2Fluid(fluidRegions.size());
PtrList<volVectorField> UFluid(fluidRegions.size());
PtrList<volScalarField> pFluid(fluidRegions.size());
PtrList<volScalarField> TFluid(fluidRegions.size());

PtrList<surfaceScalarField> phiFluid(fluidRegions.size());
PtrList<surfaceScalarField> alphaPhiFluid(fluidRegions.size());
PtrList<surfaceScalarField> rhoPhiFluid(fluidRegions.size());
PtrList<surfaceScalarField> rhoCpPhiFluid(fluidRegions.size());

PtrList<volScalarField> rhoFluid(fluidRegions.size());
PtrList<volScalarField> rhoCpFluid(fluidRegions.size());
PtrList<volScalarField> muFluid(fluidRegions.size());
PtrList<volScalarField> muMetalFluid(fluidRegions.size());
PtrList<volScalarField> kappaFluid(fluidRegions.size());
PtrList<volScalarField> epsDotFluid(fluidRegions.size());
PtrList<volScalarField> QviscFluid(fluidRegions.size());

PtrList<fv::options> fluidFvOptions(fluidRegions.size());

List<label> pRefCellFluid(fluidRegions.size(), -1);
List<scalar> pRefValueFluid(fluidRegions.size(), Zero);

forAll(fluidRegions, i)
{
    fvMesh& mesh = fluidRegions[i];

    Info<< "*** Reading fluid region " << mesh.name() << nl << endl;

    phasePropsFluid.set
    (
        i,
        new IOdictionary
        (
            IOobject
            (
                "phaseProperties",
                runTime.constant(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        )
    );

    toolPropsFluid.set
    (
        i,
        new IOdictionary
        (
            IOobject
            (
                "toolProperties",
                runTime.constant(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        )
    );

    {
        const word model
        (
            phasePropsFluid[i].subDict("metal")
           .getOrDefault<word>("viscosityModel", "constant")
        );

        if (model != "constant" && model != "SellarsTegart")
        {
            FatalIOErrorInFunction(phasePropsFluid[i])
                << "metal/viscosityModel desconocido: " << model
                << ". Opciones: constant, SellarsTegart"
                << exit(FatalIOError);
        }

        Info<< "    Metal viscosity model: " << model << nl << endl;
    }

    for (const word phase : {"metal", "air"})
    {
        const scalarList kc
        (
            phasePropsFluid[i].subDict(phase).get<scalarList>("kappaCoeffs")
        );

        if (kc.size() != 4)
        {
            FatalIOErrorInFunction(phasePropsFluid[i])
                << "kappaCoeffs de " << phase
                << " debe tener 4 coeficientes (k0 k1 k2 k3)"
                << exit(FatalIOError);
        }
    }

    // --- Campos leídos

    UFluid.set
    (
        i,
        new volVectorField
        (
            IOobject
            (
                "U",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh
        )
    );

    pFluid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "p",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh
        )
    );
    mesh.setFluxRequired(pFluid[i].name());

    alpha1Fluid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "alpha.metal",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh
        )
    );

    alpha2Fluid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "alpha.air",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            scalar(1) - alpha1Fluid[i]
        )
    );

    TFluid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "T",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            mesh
        )
    );

    // --- Flujos

    phiFluid.set
    (
        i,
        new surfaceScalarField
        (
            IOobject
            (
                "phi",
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            fvc::flux(UFluid[i])
        )
    );

    alphaPhiFluid.set
    (
        i,
        new surfaceScalarField
        (
            IOobject
            (
                "alphaPhi0.metal",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            phiFluid[i]*fvc::interpolate(alpha1Fluid[i])
        )
    );

    rhoPhiFluid.set
    (
        i,
        new surfaceScalarField
        (
            IOobject
            (
                "rhoPhi",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar(dimMass/dimTime, Zero)
        )
    );

    rhoCpPhiFluid.set
    (
        i,
        new surfaceScalarField
        (
            IOobject
            (
                "rhoCpPhi",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar(dimPower/dimTemperature, Zero)
        )
    );

    // --- Propiedades de la mezcla (se calculan en updateFluidProperties.H)

    rhoFluid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "rho",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar(dimDensity, Zero)
        )
    );

    rhoCpFluid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "rhoCp",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar(dimDensity*dimSpecificHeatCapacity, Zero)
        )
    );

    muFluid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "mu",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar(dimDynamicViscosity, Zero)
        )
    );

    // Viscosidad de la fase metal (constante o Sellars-Tegart)
    muMetalFluid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "mu.metal",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar(dimDynamicViscosity, Zero)
        )
    );

    // "kappa" lo lee la condición de acoplamiento con el sustrato
    // (kappaMethod lookup)
    kappaFluid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "kappa",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar(dimPower/dimLength/dimTemperature, Zero)
        )
    );

    epsDotFluid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "epsDot",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar(dimless/dimTime, Zero)
        )
    );

    QviscFluid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "Qvisc",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar(dimPower/dimVolume, Zero)
        )
    );

    fluidFvOptions.set(i, new fv::options(mesh));

    setRefCell
    (
        pFluid[i],
        mesh.solutionDict().subDict("PIMPLE"),
        pRefCellFluid[i],
        pRefValueFluid[i]
    );

    // Propiedades y flujos iniciales (antes de crear los sólidos, que
    // necesitan kappa del fluido en la interfaz)
    {
        #include "setRegionFluidFields.H"
        #include "updateFluidProperties.H"

        rhoPhi = phi*fvc::interpolate(rho);
        rhoCpPhi = phi*fvc::interpolate(rhoCp);
    }
}
