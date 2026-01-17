# Gravitational Hertz Experiment

This script is used to calculate the metric tensor of the gravitational wave induced by a spot-like source, which means the distance from the detector to the source is much longer than the size of the source but may be akin to the wavelength, say
\[
    R \gg d \\
    R \approx \lambda
\]

We assume that the source is a rotating column with four holes inside, where the diameter of the column is \(D \approx 5 \text{m}\) while the holes has diameter of about \(d \approx 1 \text{m}\). The columns are made of carbonfiber.

According to the linearized theory, the metric tensor at the detecter will be proportional to the second-order derivative of the quardrupole tensor,
\[
    h_{ij}^{TT} = \frac{2 G}{R c^{4}} \ddot{I}_{ij}^{TT} (t_{\text{rev}})
\]
where the quadrupole tensor is defined by
\[
    I_{ij} (t) = \int_{V'} \rho \left(x'_{i} x'_{j} - \frac{1}{3} \delta_{ij} r'^{2}\right)  \, \text{d} \tau' 
\]

For the columns under this circumstance, the calculation of the quadrupole tensor could be simplified to the linear superposition of a column with the density \(\rho\) and four columns of the density \(- \rho\).

After getting the numerical value of \(h_{ij}\), we'll use the projector to get the TT mode:
\[
    \Lambda_{ij, kl} = P_{ik} P_{jl} - \frac{1}{2} P_{ij} P_{kl}
\]
where 
\[
    P = \begin{pmatrix} 
    1 & 0 & 0 \\
    0 & 1 & 0 \\
    0 & 0 & 0
    \end{pmatrix}
\]
So \(h_{ij}^{TT} = \Lambda_{ij, kl} h_{kl}\)

The script will calculate the components of the metric tensor for a rotating carbon fiber column and output the result in the TT gauge.