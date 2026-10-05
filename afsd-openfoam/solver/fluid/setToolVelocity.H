// Velocidad en la herramienta, la barra y la superficie del sustrato
// (tesis ec. 12-13), marco de la herramienta, eje de giro = eje z por el
// origen, avance Vt según +x.
//
//   Superficie del sustrato (no deslizamiento, sólido moviéndose a -Vt):
//     u_s = -Vt
//
//   Cara de la herramienta (deslizamiento parcial delta, ec. 13):
//     u_f = ( -(1-delta) omega y - delta Vt,  (1-delta) omega x,  0 )
//   Entrada de la barra (ec. 12, literal de la tesis):
//     u_i = ( -omega y - Vt,  omega x,  -V_i )
{
    const scalar toolOmega(toolDict.get<scalar>("omega"));
    const scalar toolDelta(toolDict.get<scalar>("delta"));
    const scalar feedVi(toolDict.get<scalar>("feedVelocity"));
    const scalar traverseVt(Vt.value().x());

    const wordList toolPatches(toolDict.get<wordList>("toolPatches"));
    const word inletPatch(toolDict.get<word>("inletPatch"));

    auto checkedPatch = [&](const word& name)
    {
        const label patchi = mesh.boundaryMesh().findPatchID(name);

        if (patchi < 0)
        {
            FatalErrorInFunction
                << "toolProperties: no existe el patch " << name
                << exit(FatalError);
        }

        if (!isA<fixedValueFvPatchVectorField>(U.boundaryField()[patchi]))
        {
            FatalErrorInFunction
                << "U en el patch " << name
                << " debe ser fixedValue (lo fija el solver)"
                << exit(FatalError);
        }

        return patchi;
    };

    for (const word& name : toolPatches)
    {
        const label patchi = checkedPatch(name);
        const vectorField& Cf = mesh.Cf().boundaryField()[patchi];

        vectorField Ub(Cf.size());

        forAll(Cf, facei)
        {
            const scalar x = Cf[facei].x();
            const scalar y = Cf[facei].y();

            Ub[facei] = vector
            (
               -(1 - toolDelta)*toolOmega*y - toolDelta*traverseVt,
                (1 - toolDelta)*toolOmega*x,
                0
            );
        }

        U.boundaryFieldRef()[patchi] == Ub;
    }

    // Superficie del sustrato: sigue a Vt de frameProperties (0 en dwell y
    // fill, Vt en el avance)
    for
    (
        const word& name
      : toolDict.getOrDefault<wordList>
        (
            "substratePatches",
            wordList({"build_to_substrate"})
        )
    )
    {
        const label patchi = checkedPatch(name);
        U.boundaryFieldRef()[patchi] == -Vt.value();
    }

    {
        const label patchi = checkedPatch(inletPatch);
        const vectorField& Cf = mesh.Cf().boundaryField()[patchi];

        vectorField Ub(Cf.size());

        forAll(Cf, facei)
        {
            const scalar x = Cf[facei].x();
            const scalar y = Cf[facei].y();

            Ub[facei] = vector
            (
               -toolOmega*y - traverseVt,
                toolOmega*x,
               -feedVi
            );
        }

        U.boundaryFieldRef()[patchi] == Ub;
    }
}
