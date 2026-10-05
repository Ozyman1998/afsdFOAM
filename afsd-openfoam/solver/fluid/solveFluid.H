// Bucle PIMPLE de la región fluida. Termina al llegar a nOuterCorrectors o
// cuando el residuo inicial de p baja de outerCorrectorResidual.
for (int oCorr = 0; oCorr < nOuterCorr; ++oCorr)
{
    const bool finalIter = (oCorr == nOuterCorr - 1);

    if (finalIter)
    {
        mesh.data().setFinalIteration(true);
    }

    // Valor de la iteración anterior, necesario para subrelajar p
    // (relaxationFactors/fields/p < 1)
    p.storePrevIter();

    if (oCorr == 0)
    {
        #include "alphaEqn.H"
        #include "setToolVelocity.H"
    }

    #include "updateFluidProperties.H"

    scalar pInitialResidual = GREAT;

    {
        #include "UEqn.H"

        for (int corr = 0; corr < nCorr; ++corr)
        {
            #include "pEqn.H"
        }
    }

    #include "setToolHeat.H"
    #include "TEqn.H"
    #include "updateFluidProperties.H"

    if (finalIter)
    {
        mesh.data().setFinalIteration(false);
    }

    if (!finalIter && oCorr > 0 && pInitialResidual < outerCorrectorResidual)
    {
        Info<< "PIMPLE (" << mesh.name() << "): converged in "
            << oCorr + 1 << " outer iterations" << endl;
        break;
    }
}
