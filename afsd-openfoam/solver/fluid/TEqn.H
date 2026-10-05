// Energía en la región build (tesis ec. 10):
//   rho cp DT/Dt = div(k grad T) + 3 beta_TQ mu epsDot^2
// con epsDot = sqrt(2/3 D:D), D = symm(grad U) (ec. 7).
//
// Forma conservativa menos T*(continuidad de rhoCp). El término de
// continuidad se trata con SuSp (implícito sólo cuando refuerza la diagonal)
// para que la matriz siga acotada en la interfaz, donde rhoCp salta ~2000x.
{
    #include "calcEpsDot.H"

    Qvisc = 3.0*betaTQ*mu*sqr(epsDot);

    const volScalarField contErrRhoCp
    (
        fvc::ddt(rhoCp) + fvc::div(rhoCpPhi)
    );

    fvScalarMatrix TEqn
    (
        fvm::ddt(rhoCp, T)
      + fvm::div(rhoCpPhi, T, "div(rhoCpPhi,T)")
      - fvm::SuSp(contErrRhoCp, T)
      - fvm::laplacian(kappa, T, "laplacian(kappa,T)")
     ==
        Qvisc
      + fvOptions(rhoCp, T)
    );

    TEqn.relax();

    fvOptions.constrain(TEqn);

    TEqn.solve(T.select(finalIter));

    fvOptions.correct(T);

    Info<< "Min/max T:" << gMin(T.primitiveField()) << ' '
        << gMax(T.primitiveField()) << endl;
}
