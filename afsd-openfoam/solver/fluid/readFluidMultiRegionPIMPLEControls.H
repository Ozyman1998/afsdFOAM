    // Controles PIMPLE de la región fluida: system/<region>/fvSolution
    const dictionary& pimpleDict = mesh.solutionDict().subDict("PIMPLE");

    const int nOuterCorr =
        pimpleDict.getOrDefault<int>("nOuterCorrectors", 1);

    const int nCorr =
        pimpleDict.getOrDefault<int>("nCorrectors", 2);

    const int nNonOrthCorr =
        pimpleDict.getOrDefault<int>("nNonOrthogonalCorrectors", 0);

    const bool momentumPredictor =
        pimpleDict.getOrDefault("momentumPredictor", true);

    // Tolerancia del bucle externo: residuo inicial de p (tesis 2.2.3:
    // los sólidos se resuelven cuando el fluido alcanza la tolerancia)
    const scalar outerCorrectorResidual =
        pimpleDict.getOrDefault<scalar>("outerCorrectorResidual", 0);

    // Término transpuesto del esfuerzo viscoso (ver UEqn.H):
    //   off | mixture (grad(mu) de la mezcla) | metal (alpha*grad(mu.metal))
    const word transposeStress
    (
        pimpleDict.getOrDefault<word>("transposeStress", "off")
    );
