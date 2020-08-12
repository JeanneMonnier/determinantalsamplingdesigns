\name{determinantalsamplingdesigns}
\alias{determinantalsamplingdesigns}
\title{determinantalsamplingdesigns}

\keyword{sampling}
\keyword{design}
\keyword{determinantal}

\description{
__Introduction__

The goal of sampling theory is to acquire knowledge of a parameter of interest using only partial information. This is done by means of a sampling design, through which a random subset is observed, and the construction of an estimator. The properties of the
sampling design are thus of crucial importance to get “good” estimators. In practice, the following
issues are fundamental:

- simplicity of the design (in terms of its definition, theory and/or drawing algorithm),

- knowledge of the first and, possibly, second order inclusion probabilities,

- control of the size of the sample,

- effective construction, in particular with prescribed unequal probabilities,

- statistical amenability (consistency, central limit theorem,...),

- low Mean Square Error (MSE)/Variance of specific estimators based on the design.


__Determinantal Sampling Designs__

Determinantal Sampling Designs are a family of sampling designs that adresses all theses issues. They are indexed by Hermitian contracting matrices, called kernel, and have known inclusion probabilities for any order, among many others properties.

__This Package__

Because they have such good properties but their construction is not so easy, it is important to develop tools to make it available. This package provides two algorithms of constructions of DSD's kernel matrices and a sampling algorithm. All the results used in this package are presented in the scientific paper \code{Determinantal Sampling Designs} written by Vincent Loonis and Xavier Mary (See References).
}


\details{
The package contains three functions:

\link{loonismary} constructs a real projection matrix, kernels of fixed size sampling designs with first order inclusion probabilities given by \code{Pi}.

\link{periodicdsd} constructs a complex matrix, kernel of a determinantal sampling designs such that first order inclusion probabilities are equal and the sampling design presents periodic behaviour.

For the last two functions, it is also possible to return a matrix \code{V} such that \code{VV*} is the kernel matrix. This matrix \code{V} is an n orthonormal basis of the range of the kernel matrix.

\link{lavancier} samples from fixed size determinantal sampling designs.}


\examples{
#loonismary(c(0.5, 0.75, 0.75, 0.2, 0.4, 0.6, 0.8))
#returns a real matrix of size 7x4 with fixed first order inclusion
#probabilities equal to the ones given in input

# 0.7071068  0.0000000 0.0000000  0.0000000
# 0.5000000 -0.7071068 0.0000000  0.0000000
# 0.5000000  0.7071068 0.0000000  0.0000000
# 0.0000000  0.0000000 0.4472136  0.0000000
# 0.0000000  0.0000000 0.6324555  0.0000000
# 0.0000000  0.0000000 0.5163978 -0.5773503
# 0.0000000  0.0000000 0.3651484  0.8164966



#loonismary(rep(120/4000, 4000))
#returns a real matrix of size 4000x120 with fixed first order inclusion
#probabilities equal to 120/4000



#loonismary(rep(120/4000, 4000), TRUE)
#returns a matrix K, of size 4000x4000, kernel of the sampling design
#corresponding





#periodicdsd(4, 1, 7)
#returns a complex matrix of size 7x4 with fixed first order inclusion
#probabilities equal to 4/7

# 0.3779645+0i  0.2356570+0.2955045i -0.0841050+0.3684881i -0.3405342+0.1639926i
# 0.3779645+0i -0.0841050+0.3684881i -0.3405342-0.1639926i  0.2356570-0.2955045i
# 0.3779645+0i -0.3405342+0.1639926i  0.2356570-0.2955045i -0.0841050+0.3684881i
# 0.3779645+0i -0.3405342-0.1639926i  0.2356570+0.2955045i -0.0841050-0.3684881i
# 0.3779645+0i -0.0841050-0.3684881i -0.3405342+0.1639926i  0.2356570+0.2955045i
# 0.3779645+0i  0.2356570-0.2955045i -0.0841050-0.3684881i -0.3405342-0.1639926i
# 0.3779645+0i  0.3779645-0.0000000i  0.3779645-0.0000000i  0.3779645-0.0000000i




#periodicdsd(15, 1, 365)
#returns a complex matrix of size 7x4 with fixed first order inclusion
#probabilities equal to 15/365



#periodicdsd(15, 2, 365)
#returns an other complex matrix of size 7x4 with fixed first order inclusion
#probabilities equal to 15/365 but with different periodic characteristics.



#periodicdsd(4, 2, 7, TRUE)
#returns a matrix K, of size 7x7, kernel of the sampling design corresponding





#lavancier(V)
#returns a vector containing the elements sampled

# 1 2 6 7



#lavancier(V, B = TRUE)
#returns a vector of size equal to the number of rows of V, containing zeros
#and ones, the ones indicating the elements sampled.

# 1 0 1 0 1 1 0
#here the sample is (1, 3, 5, 6)



#lavancier(V, 5)
#returns a dataframe containing the five samples drawn

#     Sample 1 Sample 2 Sample 3 Sample 4 Sample 5
#            1        2        1        1        2
#            3        3        3        2        3
#            6        5        6        4        5
#            7        6        7        7        7


#lavancier(V, B = TRUE, C = c(FALSE, TRUE, TRUE, TRUE))
#returns a vector of size equal to the number of rows of V and plots the
#probability weights at second, third and fourth draws

# 1 1 0 0 0 1 1

}
\references{
V., Loonis, X., Mary, (2017) "Determinantal Sampling Designs".
}
\seealso{
\link{lavancier}, \link{periodicdsd}, \link{loonismary}.
}
