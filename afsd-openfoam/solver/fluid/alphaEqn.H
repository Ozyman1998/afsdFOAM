// Fracción de volumen del metal (tesis ec. 1), VOF con MULES explícito y
// compresión de interfaz. Sin tensión superficial (tesis: Ca >> 1).
{
    const dictionary& alphaControls = mesh.solverDict(alpha1.name());

    const label nAlphaCorr
    (
        alphaControls.getOrDefault<label>("nAlphaCorr", 1)
    );
    const scalar cAlpha
    (
        alphaControls.getOrDefault<scalar>("cAlpha", 1.0)
    );

    const word alphaScheme("div(phi,alpha)");
    const word alpharScheme("div(phirb,alpha)");

    // Normal a la interfaz en las caras
    const dimensionedScalar deltaN
    (
        "deltaN",
        1e-8/cbrt(average(mesh.V()))
    );
    const surfaceVectorField gradAlphaf
    (
        fvc::interpolate(fvc::grad(alpha1))
    );
    const surfaceScalarField nHatf
    (
        (gradAlphaf/(mag(gradAlphaf) + deltaN)) & mesh.Sf()
    );

    // Flujo de compresión
    surfaceScalarField phic(mag(phi/mesh.magSf()));
    phic = min(cAlpha*phic, max(phic));
    const surfaceScalarField phir(phic*nHatf);

    for (label aCorr = 0; aCorr < nAlphaCorr; ++aCorr)
    {
        alphaPhi =
            fvc::flux(phi, alpha1, alphaScheme)
          + fvc::flux
            (
                -fvc::flux(-phir, alpha2, alpharScheme),
                alpha1,
                alpharScheme
            );

        MULES::explicitSolve
        (
            geometricOneField(),
            alpha1,
            phi,
            alphaPhi,
            zeroField(),
            zeroField(),
            oneField(),
            zeroField()
        );

        alpha2 = scalar(1) - alpha1;
    }

    // Flujos de masa y de capacidad calorífica consistentes con alpha
    rhoPhi = alphaPhi*(rho1 - rho2) + phi*rho2;
    rhoCpPhi = alphaPhi*(rho1*Cp1 - rho2*Cp2) + phi*(rho2*Cp2);

    Info<< "Metal volume fraction = "
        << alpha1.weightedAverage(mesh.Vsc()).value()
        << "  Min(alpha.metal) = " << gMin(alpha1.primitiveField())
        << "  Max(alpha.metal) = " << gMax(alpha1.primitiveField())
        << endl;
}
