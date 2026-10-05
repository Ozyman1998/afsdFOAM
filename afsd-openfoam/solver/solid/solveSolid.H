{
    // AFSD (Kincaid 2022, ec. 19): el sólido se traslada a -Vt en el marco de
    // la herramienta, así que a la conducción se añade la advección:
    //     rho cp dT/dt + rho cp (-Vt).grad(T) = div(k grad(T))
    // En forma de entalpía con flujo másico de cara phiSolid = rho (-Vt).Sf
    // (rho constante por región -> phiSolid sin divergencia).
    const surfaceScalarField phiSolid
    (
        "phiSolid",
        fvc::interpolate(betav*rho)*((-Vt) & mesh.Sf())
    );

    fvScalarMatrix hEqn
    (
        fvm::ddt(betav*rho, h)
      + fvm::div(phiSolid, h, "div(phiSolid,h)")
        - (
            thermo.isotropic()
            ? fvm::laplacian(betav*thermo.alpha(), h, "laplacian(alpha,h)")
            : fvm::laplacian(betav*taniAlpha(), h, "laplacian(alpha,h)")
        )
        ==
        fvOptions(rho, h)
    );

    hEqn.relax();

    fvOptions.constrain(hEqn);

    if (coupled)
    {
        fvMatrixAssemblyPtr->addFvMatrix(hEqn);
    }
    else
    {
        Info<< "\nSolving for solid region "<< solidRegions[i].name() << endl;

        if (finalIter)
        {
            mesh.data().setFinalIteration(true);
        }

        hEqn.solve(h.select(finalIter));

        fvOptions.correct(h);

        thermo.correct();

        Info<< "Min/max T:" << min(thermo.T()).value() << ' '
            << max(thermo.T()).value() << endl;

        if (finalIter)
        {
            mesh.data().setFinalIteration(false);
        }
    }
}
