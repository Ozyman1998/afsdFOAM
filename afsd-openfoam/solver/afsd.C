/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Basado en chtMultiRegionFoam (OpenFOAM v2412, GPL v3).

Application
    afsdFoam

Description
    Additive friction stir deposition (Kincaid 2022, cap. 2).

    Todo el dominio se resuelve en el marco de la herramienta (velocidad Vt):
      - build (fluido): VOF incompresible metal + aire, PIMPLE, energía con
        conducción k(T) y calentamiento viscoso.
      - substrate, anvil (sólidos): conducción + advección a -Vt.

    Algoritmo por paso de tiempo (tesis 2.2.3):
      1. PIMPLE en build hasta tolerancia (outerCorrectorResidual)
      2. Energía en sólidos
      3. nEnergyCouplingLoops bucles extra solo de energía (build con flujo
         congelado + sólidos)

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "turbulentFluidThermoModel.H"
#include "rhoReactionThermo.H"
#include "CombustionModel.H"
#include "fixedGradientFvPatchFields.H"
#include "fixedValueFvPatchFields.H"
#include "mixedFvPatchFields.H"
#include "regionProperties.H"
#include "solidRegionDiffNo.H"
#include "solidThermo.H"
#include "radiationModel.H"
#include "fvOptions.H"
#include "coordinateSystem.H"
#include "loopControl.H"
#include "pressureControl.H"
#include "MULES.H"
#include "geometricOneField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "AFSD: VOF incompresible metal/aire con energía, acoplado a sólidos"
        " con advección, en el marco de la herramienta."
    );

    #define NO_CONTROL
    #define CREATE_MESH createMeshesPostProcess.H
    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMeshes.H"
    #include "createFields.H"
    #include "createTimeControls.H"
    #include "readSolidTimeControls.H"
    #include "incompressibleMultiRegionCourantNo.H"
    #include "solidRegionDiffusionNo.H"
    #include "setInitialMultiRegionDeltaT.H"

    // Acoplamiento térmico secuencial entre regiones (tesis), sin ensamblado
    // implícito conjunto. solveSolid.H sigue refiriéndose a estas variables.
    const bool coupled = false;
    autoPtr<fvMatrix<scalar>> fvMatrixAssemblyPtr;

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readTimeControls.H"
        #include "readSolidTimeControls.H"

        // Controles globales (system/fvSolution, PIMPLE)
        fvSolution topSolution(runTime);
        const label nEnergyCouplingLoops =
            topSolution.subDict("PIMPLE")
           .getOrDefault<label>("nEnergyCouplingLoops", 0);

        #include "incompressibleMultiRegionCourantNo.H"
        #include "solidRegionDiffusionNo.H"
        #include "setMultiRegionDeltaT.H"

        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- 1. Regiones fluidas: PIMPLE hasta tolerancia
        forAll(fluidRegions, i)
        {
            fvMesh& mesh = fluidRegions[i];

            Info<< "\nSolving for fluid region " << mesh.name() << endl;

            #include "setRegionFluidFields.H"
            #include "readFluidMultiRegionPIMPLEControls.H"
            #include "solveFluid.H"
        }

        // --- 2. Regiones sólidas
        {
            const bool finalIter = (nEnergyCouplingLoops == 0);

            forAll(solidRegions, i)
            {
                fvMesh& mesh = solidRegions[i];

                #include "readSolidMultiRegionPIMPLEControls.H"
                #include "setRegionSolidFields.H"
                #include "solveSolid.H"
            }
        }

        // --- 3. Bucles extra solo de energía
        for (label eLoop = 0; eLoop < nEnergyCouplingLoops; ++eLoop)
        {
            const bool finalIter = (eLoop == nEnergyCouplingLoops - 1);

            Info<< "\nEnergy coupling loop " << eLoop + 1 << " of "
                << nEnergyCouplingLoops << endl;

            forAll(fluidRegions, i)
            {
                fvMesh& mesh = fluidRegions[i];

                Info<< "\nSolving energy for fluid region " << mesh.name()
                    << endl;

                #include "setRegionFluidFields.H"
                #include "setToolHeat.H"
                #include "TEqn.H"
                #include "updateFluidProperties.H"
            }

            forAll(solidRegions, i)
            {
                fvMesh& mesh = solidRegions[i];

                #include "readSolidMultiRegionPIMPLEControls.H"
                #include "setRegionSolidFields.H"
                #include "solveSolid.H"
            }
        }

        runTime.write();

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
