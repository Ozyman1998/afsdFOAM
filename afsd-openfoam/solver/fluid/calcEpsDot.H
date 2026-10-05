// Tasa de deformación efectiva (tesis ec. 7): epsDot = sqrt(2/3 D:D),
// D = symm(grad U)
{
    const volSymmTensorField D(symm(fvc::grad(U)));
    epsDot = sqrt((2.0/3.0)*(D && D));
}
