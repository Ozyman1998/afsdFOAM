/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  v2412                                 |
|   \\  /    A nd           |                                                 |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       volScalarField;
    location    "0/build";
    object      alpha.metal;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Fracción de volumen de metal (1 = AA6061, 0 = aire)

dimensions      [0 0 0 0 0 0 0];

internalField   uniform 0;

boundaryField
{
    "build_(top|sides|front|back)"
    {
        type            inletOutlet;
        inletValue      uniform 0;
        value           uniform 0;
    }

    "build_to_.*"
    {
        type            zeroGradient;
    }

    "tool(Face|Side)"
    {
        type            zeroGradient;
    }

    inlet
    {
        type            fixedValue;
        value           uniform 1;
    }
}

// ************************************************************************* //
