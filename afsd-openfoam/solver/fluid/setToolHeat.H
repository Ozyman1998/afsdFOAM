// Condiciones térmicas de la herramienta y de la barra.
//
// Cara de la herramienta (tesis ec. 15-17):
//   q_i     = [ (1-delta) beta_TQ tau + delta mu_f p ] alpha omega r
//   k dT/dn = f q_i,     f = J_B/(J_B+J_T) (fracción hacia el metal)
//   tau     = sigma_f(T, epsDot = tauStrainRate)   (Sellars-Tegart, ec. 8)
//   mu_f    = muF0 exp(-lambda delta omega r)      (Nandan; cap. 3)
//   p       = presión local en la cara, limitada a p >= 0
//
// Entrada de la barra (ec. 14): conducción 1D hasta Ta a una distancia L
//   k dT/dn = -(k/L)(T - Ta)  ->  dT/dn = -(T - Ta)/L   (k se cancela)
// implementado como condición mixta con valueFraction = 1/(1 + L*deltaCoeffs)
{
    if (!isSellarsTegart)
    {
        FatalErrorInFunction
            << "El calor de la herramienta usa tau = sigma_f(T) de "
            << "Sellars-Tegart: metal/viscosityModel debe ser SellarsTegart"
            << exit(FatalError);
    }

    const scalar toolOmega(toolDict.get<scalar>("omega"));
    const scalar toolDelta(toolDict.get<scalar>("delta"));
    const scalar muF0(toolDict.get<scalar>("muF0"));
    const scalar muFLambda(toolDict.get<scalar>("lambda"));
    const scalar heatFraction(toolDict.get<scalar>("heatFraction"));
    const scalar tauEpsDot(toolDict.get<scalar>("tauStrainRate"));
    const scalar rodL(toolDict.get<scalar>("rodLength"));
    const scalar rodTa(toolDict.get<scalar>("rodTa"));

    const wordList toolPatches(toolDict.get<wordList>("toolPatches"));
    const word inletPatch(toolDict.get<word>("inletPatch"));

    // Tensión de fluencia de Sellars-Tegart a la tasa tauStrainRate
    auto flowStress = [&](const scalar Tval)
    {
        const scalar Te = min(max(Tval, stTmin), stTmax);
        const scalar Z = tauEpsDot*std::exp(stQ/(stR*Te));
        return stSigmaN*std::asinh(std::pow(Z/stA, 1.0/stN));
    };

    scalar toolHeat = 0;   // [W] hacia el metal
    scalar rodHeat = 0;    // [W] que sale por la barra

    for (const word& name : toolPatches)
    {
        const label patchi = mesh.boundaryMesh().findPatchID(name);

        if
        (
            patchi < 0
         || !isA<fixedGradientFvPatchScalarField>(T.boundaryField()[patchi])
        )
        {
            FatalErrorInFunction
                << "T en el patch " << name
                << " debe existir y ser fixedGradient (lo fija el solver)"
                << exit(FatalError);
        }

        fixedGradientFvPatchScalarField& Tp =
            refCast<fixedGradientFvPatchScalarField>
            (
                T.boundaryFieldRef()[patchi]
            );

        const vectorField& Cf = mesh.Cf().boundaryField()[patchi];
        const scalarField& magSf = mesh.magSf().boundaryField()[patchi];
        const scalarField& pb = p.boundaryField()[patchi];
        const scalarField& ab = alpha1.boundaryField()[patchi];
        const scalarField& kb = kappa.boundaryField()[patchi];

        scalarField& gradT = Tp.gradient();

        forAll(gradT, facei)
        {
            const scalar r = Foam::sqrt(sqr(Cf[facei].x()) + sqr(Cf[facei].y()));
            const scalar a = min(max(ab[facei], scalar(0)), scalar(1));
            const scalar tau = flowStress(Tp[facei]);
            const scalar muF = muF0*std::exp(-muFLambda*toolDelta*toolOmega*r);

            const scalar q =
                (
                    (1 - toolDelta)*betaTQ*tau
                  + toolDelta*muF*max(pb[facei], scalar(0))
                )*a*toolOmega*r;

            gradT[facei] = heatFraction*q/max(kb[facei], SMALL);
            toolHeat += heatFraction*q*magSf[facei];
        }
    }

    {
        const label patchi = mesh.boundaryMesh().findPatchID(inletPatch);

        if
        (
            patchi < 0
         || !isA<mixedFvPatchScalarField>(T.boundaryField()[patchi])
        )
        {
            FatalErrorInFunction
                << "T en el patch " << inletPatch
                << " debe existir y ser mixed (lo fija el solver)"
                << exit(FatalError);
        }

        mixedFvPatchScalarField& Ti =
            refCast<mixedFvPatchScalarField>(T.boundaryFieldRef()[patchi]);

        Ti.refValue() = rodTa;
        Ti.refGrad() = 0;
        Ti.valueFraction() = 1.0/(1.0 + rodL*Ti.patch().deltaCoeffs());

        const scalarField& magSf = mesh.magSf().boundaryField()[patchi];
        const scalarField& kb = kappa.boundaryField()[patchi];

        forAll(Ti, facei)
        {
            rodHeat += kb[facei]*(Ti[facei] - rodTa)/rodL*magSf[facei];
        }
    }

    reduce(toolHeat, sumOp<scalar>());
    reduce(rodHeat, sumOp<scalar>());

    Info<< "Tool heat into build [W] = " << toolHeat
        << "   Heat lost up the rod [W] = " << rodHeat << endl;
}
