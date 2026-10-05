// Propiedades de la mezcla: media ponderada por volumen (tesis ec. 2).
// rhoCp se promedia como producto (rho*Cp) de cada fase para que la ecuación
// de energía sea conservativa con el flujo rhoCpPhi.
{
    const volScalarField a1
    (
        min
        (
            max(alpha1, dimensionedScalar(dimless, Zero)),
            dimensionedScalar(dimless, 1.0)
        )
    );
    const volScalarField a2(scalar(1) - a1);

    rho == a1*rho1 + a2*rho2;
    rhoCp == a1*(rho1*Cp1) + a2*(rho2*Cp2);

    // --- Viscosidad del metal
    if (isSellarsTegart)
    {
        #include "calcEpsDot.H"

        // Sellars-Tegart / Sheppard-Wright (tesis ec. 6, 8, 9):
        //   Z     = epsDot exp(Q/(R T))
        //   sigma = sigmaN asinh( (Z/A)^(1/n) )
        //   mu    = sigma/(3 epsDot)
        // Límites numéricos: epsDot >= epsDotMin, Tmin <= T <= Tmax,
        // mu <= muMax.
        auto sellarsTegartMu = [&](const scalar Tval, const scalar epsVal)
        {
            const scalar e = max(epsVal, stEpsDotMin);
            const scalar Te = min(max(Tval, stTmin), stTmax);
            const scalar Z = e*std::exp(stQ/(stR*Te));
            const scalar sigma =
                stSigmaN*std::asinh(std::pow(Z/stA, 1.0/stN));
            return min(sigma/(3.0*e), stMuMax);
        };

        {
            scalarField& muI = muMetal.primitiveFieldRef();
            const scalarField& TI = T.primitiveField();
            const scalarField& eI = epsDot.primitiveField();

            forAll(muI, celli)
            {
                muI[celli] = sellarsTegartMu(TI[celli], eI[celli]);
            }
        }

        forAll(muMetal.boundaryField(), patchi)
        {
            scalarField& muP = muMetal.boundaryFieldRef()[patchi];
            const scalarField& TP = T.boundaryField()[patchi];
            const scalarField& eP = epsDot.boundaryField()[patchi];

            forAll(muP, facei)
            {
                muP[facei] = sellarsTegartMu(TP[facei], eP[facei]);
            }
        }
    }
    else
    {
        muMetal == mu1;
    }

    mu == a1*muMetal + a2*mu2;

    // --- k = k0 + k1 T + k2 T^2 + k3 T^3, T en °C (tesis ec. 11).
    // La cúbica sólo vale en su rango de ajuste (fuera puede dar k < 0):
    // se evalúa con T limitada a [kappaTmin, kappaTmax] °C.
    const volScalarField Tc
    (
        min
        (
            max
            (
                (T - dimensionedScalar(dimTemperature, 273.15))
               /dimensionedScalar(dimTemperature, 1.0),
                dimensionedScalar(dimless, kappaTmin)
            ),
            dimensionedScalar(dimless, kappaTmax)
        )
    );

    const dimensionedScalar kUnit(dimPower/dimLength/dimTemperature, 1.0);

    kappa ==
        kUnit
       *(
            a1
           *(
               kappaCoeffs1[0]
             + Tc*(kappaCoeffs1[1] + Tc*(kappaCoeffs1[2] + Tc*kappaCoeffs1[3]))
            )
          + a2
           *(
               kappaCoeffs2[0]
             + Tc*(kappaCoeffs2[1] + Tc*(kappaCoeffs2[2] + Tc*kappaCoeffs2[3]))
            )
        );
}
