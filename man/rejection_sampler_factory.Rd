% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rejections_sampler.R
\name{rejection_sampler_factory}
\alias{rejection_sampler_factory}
\title{Compile rejection samplers for arbitrary envelope}
\usage{
rejection_sampler_factory(enve, evalmode = 2)
}
\arguments{
\item{enve}{An object of class 'Envelope'.}

\item{alpha}{A scalar value ensuring \eqn{\alpha f \leq g}.}
}
\value{
A rejection sampler for the given envelope.
}
\description{
Compile rejection samplers for arbitrary envelope
}
\examples{
#Normal Distribution
enve <- LogLinearEnvelope(get_rv(), c(-2,0,1))
sampler <- rejection_sampler_factory(enve)
new_sampler <- rejection_sampler_factory(enve)
sim <- sampler(50000)
plot(enve)
sampler(5000000, adapt_enve = T)
sampler \%>\% environment \%>\% `[[`("enve") \%>\% plot()
hist(sim, prob=TRUE)
curve(enve$base_rv$f(x), -3, 3, col="blue", add=TRUE)
qqnorm(sim)
qqline(sim)
shapiro.test(sim[1:5000])
environment(sampler)$p
environment(sampler)$credibility

#Gamma Distribution
enve <- LogLinearEnvelope(get_rv("g"))
sampler <- rejection_sampler_factory(enve)
sim <- sampler(50000)
hist(sim, prob=TRUE, ylim = c(0,0.5))
curve(enve$base_rv$f(x), 0, 8, col="blue", add=TRUE)

# Laplacian Envelope
enve <- LaplaceEnvelope(dnorm, sim_method = 2)
sampler <- rejection_sampler_factory(enve)
sim <- sampler(50000)
hist(sim, prob=TRUE, ylim=c(0,0.6))
curve(enve$base_rv$f(x), -3, 3, col="blue", add=TRUE)
}
