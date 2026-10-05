// Corrector de presión (tesis ec. 3: div(U) = 0)
{
    volScalarField rAU("rAU", 1.0/UEqn.A());
    surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU));
    volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));

    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        fvc::flux(HbyA)
      + fvc::interpolate(rho*rAU)*fvc::ddtCorr(U, phi)
    );

    adjustPhi(phiHbyA, U, p);

    for (int nonOrth = 0; nonOrth <= nNonOrthCorr; ++nonOrth)
    {
        fvScalarMatrix pEqn
        (
            fvm::laplacian(rAUf, p) == fvc::div(phiHbyA)
        );

        pEqn.setReference(pRefCell, getRefCellValue(p, pRefCell));

        const scalar initialResidual = pEqn.solve
        (
            p.select
            (
                oCorr == nOuterCorr - 1
             && corr == nCorr - 1
             && nonOrth == nNonOrthCorr
            )
        ).initialResidual();

        if (corr == 0 && nonOrth == 0)
        {
            pInitialResidual = initialResidual;
        }

        if (nonOrth == nNonOrthCorr)
        {
            phi = phiHbyA - pEqn.flux();

            p.relax();

            U = HbyA - rAU*fvc::reconstruct(pEqn.flux()/rAUf);
            U.correctBoundaryConditions();
            fvOptions.correct(U);
        }
    }

    if (p.needReference())
    {
        p += dimensionedScalar
        (
            "p",
            p.dimensions(),
            pRefValue - getRefCellValue(p, pRefCell)
        );
    }

    // Errores de continuidad
    {
        const volScalarField contErr(fvc::div(phi));

        const scalar sumLocalContErr =
            runTime.deltaTValue()
           *mag(contErr)().weightedAverage(mesh.V()).value();

        Info<< "time step continuity errors (" << mesh.name() << ")"
            << ": sum local = " << sumLocalContErr << endl;
    }
}
