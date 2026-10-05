// Cantidad de movimiento (tesis ec. 4), sin gravedad ni tensión superficial:
//   rho DU/Dt = -grad(p) + div( mu (grad(U) + grad(U)^T) )
//
// Término transpuesto: con div(U) = 0,
//   div(mu grad(U)^T) = grad(U)^T . grad(mu)   [= (grad(U) & grad(mu))]
// Es explícito y sólo actúa donde mu varía. PIMPLE/transposeStress:
//   off     : no se incluye
//   mixture : con mu de la mezcla (inestable en la interfaz, paso 4)
//   metal   : alpha*grad(U)^T.grad(mu.metal), sin el salto de la interfaz
fvVectorMatrix UEqn
(
    fvm::ddt(rho, U)
  + fvm::div(rhoPhi, U, "div(rhoPhi,U)")
  - fvm::laplacian(mu, U, "laplacian(mu,U)")
 ==
    fvOptions(rho, U)
);

if (transposeStress == "mixture" || transposeStress == "on")
{
    UEqn -= (fvc::grad(U) & fvc::grad(mu));
}
else if (transposeStress == "metal")
{
    UEqn -= alpha1*(fvc::grad(U) & fvc::grad(muMetal));
}
else if (transposeStress != "off" && transposeStress != "none")
{
    FatalErrorInFunction
        << "PIMPLE/transposeStress desconocido: " << transposeStress
        << ". Opciones: off, mixture, metal" << exit(FatalError);
}

UEqn.relax();

fvOptions.constrain(UEqn);

if (momentumPredictor)
{
    solve
    (
        UEqn
     ==
        fvc::reconstruct(-fvc::snGrad(p)*mesh.magSf())
    );

    fvOptions.correct(U);
}
