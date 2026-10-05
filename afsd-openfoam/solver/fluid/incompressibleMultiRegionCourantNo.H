// Número de Courant de las regiones fluidas (flujo volumétrico, incompresible)
scalar CoNum = -GREAT;

forAll(fluidRegions, regioni)
{
    const fvMesh& regionMesh = fluidRegions[regioni];

    const scalarField sumPhi
    (
        fvc::surfaceSum(mag(phiFluid[regioni]))().primitiveField()
    );

    const scalar regionCoNum =
        0.5*gMax(sumPhi/regionMesh.V().field())*runTime.deltaTValue();

    const scalar meanCoNum =
        0.5*(gSum(sumPhi)/gSum(regionMesh.V().field()))*runTime.deltaTValue();

    Info<< "Region: " << regionMesh.name()
        << " Courant Number mean: " << meanCoNum
        << " max: " << regionCoNum << endl;

    CoNum = max(CoNum, regionCoNum);
}
